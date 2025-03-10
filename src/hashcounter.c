#include "hashcounter.h"
#include "kseq.h"

/*
 * Implemented hashtable has open addressing with linear probe collision policy
 * A memory efficient hashtable entry is also available, with half the memory usage (collisions are more likely with this)
 * The hashtable resizes when the load factor exceeds 0.75 after entering the k-mers of a subcontig
 * Keys are k-mers hashed using the non-cryptographic MurMurHash
 * Values are the status of the k-mer (i.e. unique or not) and also the id of the subcontig from which it originates
 */

KSEQ_INIT(gzFile, gzread)

hashtable_t* hashtable_create(uint32_t kmer_size, uint32_t num_subconts, uint32_t num_threads){
    hashtable_t* ht = (hashtable_t*) malloc(sizeof(hashtable_t));
    ht->subcontig_names = calloc(num_subconts,sizeof(char*));
    ht->subcontig_counts = calloc(num_subconts,sizeof(int));
    ht->num_subcontigs = num_subconts;
    ht->curr_subcontig = 1;
    ht->size = INITIAL_HT_SIZE;
    ht->count = 0;
    ht->entry_bitmask = INITIAL_HT_BITMASK;
    ht->kmer_size = kmer_size;
    ht->num_threads = num_threads;
    ht->items = (ht_element_t*) calloc(INITIAL_HT_SIZE, sizeof(ht_element_t));
    if(ht->items == NULL){
        exit(0);
    }
    sem_init(&(ht->can_edit_ht), 0, num_threads);
    pthread_mutex_init(&ht->resize_lock, NULL);
    return ht;
}

void hashtable_destroy(hashtable_t* ht){
    int i=1;
    while(ht->subcontig_names[i]!=NULL){
        free(ht->subcontig_names[i]);
        ++i;
    }
    free(ht->subcontig_names);
    free(ht->subcontig_counts);
    free(ht->items);
    sem_destroy(&(ht->can_edit_ht));
    pthread_mutex_destroy(&ht->resize_lock);
    free(ht);
}

// insert into hashtable with linear probe collision policy
// return entry if found in ht 
ht_element_t* hashtable_insert(hashtable_t* ht, uint64_t key, ht_element_status status, uint32_t subcontig_id){
    uint64_t hash = key & ht->entry_bitmask;
    ht_element_t* current_item = &(ht->items[hash]);
    ht_element_status emp = EMPTY;
    while(!atomic_compare_exchange_strong(&(current_item->status), &emp, status)){
        emp = EMPTY;
        // ensure key isn't 0 to prevent race condition where item is actively being created by another thread
        volatile uint64_t* curr_key = &(current_item->key);
        while(*curr_key == 0){}
        if(*curr_key == key) return current_item;
        ++current_item;
        // reset to beginning of hashtable if end is reached
        if((uint64_t)(current_item - ht->items) == ht->size) current_item = ht->items;
    }
    // important to write key after id, as it is assumed a filled in key means a filled in id
    current_item->subcontig_id = subcontig_id;
    current_item->key = key;
    return NULL;
}

// double ht size and re-enter all elements from left to right
void hashtable_resize(hashtable_t* ht){
    //uint64_t putative_count = 0;
    ht->size *= 2;
    printf("Hashtable is resizing, new size will use ~ %ld GiB of memory\n", ht->size/INITIAL_HT_SIZE/2);
    uint64_t changed_bit = ht->entry_bitmask;
    ht->entry_bitmask = (ht->entry_bitmask << 1) | 0x1;
    changed_bit ^= ht->entry_bitmask;
    ht->items = realloc(ht->items, ht->size * sizeof(ht_element_t));
    if (!ht->items) {
        fprintf(stderr, "Error: not enough space for hashtable to resize. Program exiting.\n");
        exit(EXIT_FAILURE);
    }
    memset(&ht->items[ht->size/2], 0, ht->size/2);
    ht_element_t* current_entry =  ht->items-1;
    while(current_entry != &ht->items[ht->size/2]){
        ++current_entry;
        if(current_entry->status == EMPTY) continue;
        if(hashtable_insert(ht, current_entry->key, current_entry->status, current_entry->subcontig_id)!=NULL) continue;
        current_entry->key = 0;
        current_entry->status = EMPTY;
        current_entry->subcontig_id = 0;
    }
}

// return the sum of all unique hashes
uint64_t sum_unique_hahses(hashtable_t* ht){
    uint64_t sum = 0;
    for(uint32_t i=0; i<ht->num_subcontigs; ++i){
        sum+=ht->subcontig_counts[i];
    }
    return sum;
}

// add one k-mer (to be hashed) to the hashtable
static inline void hashtable_add_kmer(hashtable_t* ht, char* seq, uint32_t subcont_id, int32_t* uniq_count, uint64_t* count){
    uint64_t hash = MurmurHash64A(seq, ht->kmer_size, (uint64_t)07062024);
    ht_element_t* hashtable_item = hashtable_insert(ht, hash, UNIQUE, subcont_id);
    ht_element_status uniq = UNIQUE;
    if(hashtable_item == NULL){
        ++uniq_count[subcont_id];
        ++*count;
    } else if(atomic_compare_exchange_strong(&(hashtable_item->status), &uniq, NON_UNIQUE)){
        --uniq_count[hashtable_item->subcontig_id];
    }
}


// function for marking a k-mer as non-unique
static inline void hashtable_mark_kmer(hashtable_t* ht, char* seq, uint32_t subcont_id, int32_t* uniq_count, uint64_t* count){
    uint64_t hash = MurmurHash64A(seq, ht->kmer_size, (uint64_t)07062024);
    if(hashtable_insert(ht, hash, NON_UNIQUE, subcont_id)==NULL) ++*count;
}

/* following function adapted from Austin Appleby */
uint64_t MurmurHash64A (const void* key, int len, uint64_t seed){
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;
    uint64_t h = seed ^ (len * m);
    const uint64_t* data = (const uint64_t*) key;
    const uint64_t* end = data + (len/8);

    while(data != end){
        uint64_t k = *data++;

        k *= m; 
        k ^= k >> r; 
        k *= m; 

        h ^= k;
        h *= m; 
    }
    const unsigned char * data2 = (const unsigned char*)data;
    switch(len & 7){
        case 7: h ^= (uint64_t)(data2[6]) << 48; break;
        case 6: h ^= (uint64_t)(data2[5]) << 40; break;
        case 5: h ^= (uint64_t)(data2[4]) << 32; break;
        case 4: h ^= (uint64_t)(data2[3]) << 24; break;
        case 3: h ^= (uint64_t)(data2[2]) << 16; break;
        case 2: h ^= (uint64_t)(data2[1]) << 8; break;
        case 1: h ^= (uint64_t)(data2[0]); break;
            h *= m;
    };
    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

// following lookup basemap for use in reverse complementing
static const unsigned char basemap[256] = {
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     96, 'T', 'v', 'G', 'h', 'e', 'f', 'C', 'd', 'i', 'j', 'm', 'l', 'k', 'N', 'o',
    'p', 'q', 'y', 's', 'A', 'A', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};
static inline char* reverse_complement(char* seq){
    char* rc = calloc(strlen(seq)+1, sizeof(char));
    uint32_t j = 0;
    for (int i = strlen(seq) - 1; i >= 0; --i) {
        rc[j] = basemap[(int)seq[i]];
        ++j;
    }
    return rc;
}

// check if sequence k-mer has an N in it
static inline int check_n(char* seq, uint32_t kmer_size){
    for(uint32_t i=0; i < kmer_size; ++i){
        if(seq[i] == 'N') return i+1;
    }
    return 0;
}

// return kseq to subcontig sequence given directory and file name
// caller is responsible for destroying kseq
kseq_t* open_subcontig(hashtable_t* ht, char* dir_location, char* subcont_name, uint32_t subcont_id){
    kseq_t* seq;
    char* subcont_location;
    uint32_t loc_size = strlen(dir_location)+strlen(subcont_name)+1;
    subcont_location = calloc(loc_size, sizeof(char));
    strcpy(subcont_location, dir_location);
    subcont_location = strcat(subcont_location, subcont_name);
    gzFile fp = gzopen(subcont_location,"r");
    if(fp == NULL){
        fprintf(stderr, "Error opening %s\n", subcont_name);
        exit(EXIT_FAILURE);
    }
    seq = kseq_init(fp);
    kseq_read(seq);
    gzclose(fp);
    return(seq);
}

// threads start here via pthread_create
void* start_hash_thread(void* input){
    job_queue_t* job_queue = (job_queue_t*) input;
    hashtable_t* original_ht = job_queue->job_head->ht;
    if(!original_ht) return NULL;
    // counts shall be specific to this thread until the end where they are atomically combined
    int32_t* thread_subcontig_counts = calloc(job_queue->job_head->ht->num_subcontigs, sizeof(int32_t));
    // job queue is shared between all threads and must be thread safe
    pthread_mutex_lock(&job_queue->job_queue_lock);
    while(job_queue->job_head->ht != NULL){
        hash_job_t* curr_job = job_queue->job_head;
        curr_job->local_uniq_counts = thread_subcontig_counts;
        ++job_queue->job_head;
        pthread_mutex_unlock(&job_queue->job_queue_lock);
        hash_and_insert_subcontig(curr_job);
        pthread_mutex_lock(&job_queue->job_queue_lock);
    }
    pthread_mutex_unlock(&job_queue->job_queue_lock);
    // update unique counts
    _Atomic int32_t* uniq_counts = original_ht->subcontig_counts;
    for(int i = 0; i < original_ht->num_subcontigs; ++i){
        uniq_counts[i] += thread_subcontig_counts[i];
    }
    free(thread_subcontig_counts);
    return NULL;
}

// add k-mers to the hashtable for an entire subcontig
void hash_and_insert_subcontig(hash_job_t* hash_job){
    uint32_t i = 0;
    hashtable_t* ht = hash_job->ht;
    kseq_t* kseq = open_subcontig(hash_job->ht, hash_job->dir_location, hash_job->subcontig_loc, hash_job->subcontig_id);
    char* subcont_header;
    // comments are technically anything after a semicolon, but will be included as part of subcont_header
    if(kseq->comment.s != NULL){
        subcont_header = calloc(strlen(kseq->name.s)+strlen(kseq->comment.s)+2, sizeof(char));
        memcpy(subcont_header, kseq->name.s, strlen(kseq->name.s));
        subcont_header[strlen(kseq->name.s)] = ' ';
        strcat(subcont_header, kseq->comment.s);
    } else {
        subcont_header = calloc(strlen(kseq->name.s)+1, sizeof(char));
        memcpy(subcont_header, kseq->name.s, strlen(kseq->name.s));
    }
    ht->subcontig_names[hash_job->subcontig_id] = subcont_header;

    char* seq = kseq->seq.s;
    char* rc = reverse_complement(seq);
    uint32_t n;
    uint32_t seq_len = strlen(seq);
    uint64_t local_counts = 0;
    while((n=check_n(&seq[i], ht->kmer_size))){
        i+=n;
    }
    sem_wait(&(ht->can_edit_ht));
    while(ht->kmer_size + i <= seq_len){
        if(seq[i+ht->kmer_size-1]=='N'){
            i+=ht->kmer_size;
            if(i+ht->kmer_size > seq_len) break;
            while((n=check_n(&seq[i], ht->kmer_size))) i+=n;
            continue;
        }
    
        // only add canonical k-mers
        if(strncmp(&seq[i], &rc[seq_len-ht->kmer_size-i], ht->kmer_size) < 0){
            hash_job->kmer_func(ht, &seq[i], hash_job->subcontig_id, hash_job->local_uniq_counts, &local_counts);
        }else{
            hash_job->kmer_func(ht, &rc[seq_len-ht->kmer_size-i], hash_job->subcontig_id, hash_job->local_uniq_counts, &local_counts);
        }
        ++i;
    }
    ht->count += local_counts;
    sem_post(&(ht->can_edit_ht));
    kseq_destroy(kseq);
    free(rc);

    // resize hashtable if load factor is >0.75 after subcontig addition
    pthread_mutex_lock(&ht->resize_lock);
    if((float) ht->count / ht->size > 0.75){
        for(int j = 0; j < ht->num_threads; ++j){
            sem_wait(&(ht->can_edit_ht));
        }
        /*
         * could multithread the below func
         */
        hashtable_resize(ht);

        for(int j = 0; j < ht->num_threads; ++j){
            sem_post(&(ht->can_edit_ht));
        }
    }
    pthread_mutex_unlock(&ht->resize_lock);
}

// add k-mers to the hashtable for all subcontigs in a directory
void hash_and_insert(hashtable_t* ht, char* dir_location, kmer_func_t kmer_func){
    struct dirent *de;
    DIR *dr = opendir(dir_location);
    job_queue_t* job_queue = malloc(sizeof(job_queue_t));
    job_queue->job_head = calloc(ht->num_subcontigs+1, sizeof(hash_job_t));
    pthread_mutex_init(&job_queue->job_queue_lock, NULL);
    uint32_t queue_start = ht->curr_subcontig;
    hash_job_t* original_job_head = job_queue->job_head;

    if(dr == NULL) {
        fprintf(stderr, "Could not open subcontigs directory\n\n");
        exit(EXIT_FAILURE);
    }
    while (((de = readdir(dr)) != NULL)) {
        if(!(strlen(de->d_name) >= 10 && strcmp(&de->d_name[strlen(de->d_name) - 10], ".subcontig") == 0)) continue;
        job_queue->job_head[ht->curr_subcontig].ht = ht;
        job_queue->job_head[ht->curr_subcontig].dir_location = dir_location;
        job_queue->job_head[ht->curr_subcontig].subcontig_loc = calloc(strlen(de->d_name)+1, sizeof(char));
        memcpy(job_queue->job_head[ht->curr_subcontig].subcontig_loc, de->d_name, strlen(de->d_name));
        job_queue->job_head[ht->curr_subcontig].subcontig_id = ht->curr_subcontig;
        job_queue->job_head[ht->curr_subcontig].kmer_func = kmer_func;
        ++ht->curr_subcontig;
    }
    
    job_queue->job_head += queue_start;
    // start threads
    pthread_t* tids = calloc(ht->num_threads, sizeof(pthread_t));
    pthread_mutex_lock(&job_queue->job_queue_lock);
    for(int i = 0; i < ht->num_threads; ++i){
        pthread_create(&tids[i], NULL, start_hash_thread, job_queue);
    }
    pthread_mutex_unlock(&job_queue->job_queue_lock);

    // wait for all threads to finish
    for(int i = 0; i < ht->num_threads; ++i){
        pthread_join(tids[i], NULL);
    }

    pthread_mutex_destroy(&job_queue->job_queue_lock);
    free(original_job_head);
    free(job_queue);
    free(tids);
    closedir(dr);
}

int main(int argc, char **argv){
    int opt;
    char* subcontigs = NULL;
    char* exc_subcontigs = NULL;
    char* outdir = NULL;
    uint32_t kmer_size = 0;
    uint32_t num_subcontigs = 0;
    uint32_t num_threads = 8;

    // parse options
    while ((opt = getopt(argc, argv, "s:e:k:n:o:t:h")) != -1) {
        switch (opt) {
            case 's': {
                subcontigs = calloc(strlen(optarg) + 2, sizeof(char));
                strcpy(subcontigs, optarg);
                subcontigs[strlen(optarg)] = '/';
            } break;
            case 'e': {
                exc_subcontigs = calloc(strlen(optarg) + 2, sizeof(char));
                strcpy(exc_subcontigs, optarg);
                exc_subcontigs[strlen(optarg)] = '/';
            } break;
            case 'k': {
                kmer_size = atoi(optarg);
            } break;
            case 'n': {
                num_subcontigs = atoi(optarg);
            } break;
            case 't': {
                num_threads = atoi(optarg);
            } break;
            case 'o': {
                outdir = calloc(strlen(optarg) + strlen("/KmerContent.report") + 1, sizeof(char));
                strcpy(outdir, optarg);
                strcat(outdir, "/KmerContent.report");
            } break;
            case 'h': {
                printf(USAGE);
                return EXIT_SUCCESS;
            }
            default: {
                printf(USAGE);
                return EXIT_FAILURE;
            }
        }
    }

    // check validity of inputs
    if(subcontigs == NULL || exc_subcontigs == NULL || outdir == NULL || kmer_size == 0 || num_subcontigs == 0) {
        printf(USAGE);
        return EXIT_FAILURE;
    }

    printf("Hashing and counting k-mers\n");
    hashtable_t* ht = hashtable_create(kmer_size, num_subcontigs+1, num_threads);

    // main pipeline
    
    printf("Hashing excluded subcontigs and marking them as non-unique\n");
    hash_and_insert(ht, exc_subcontigs, hashtable_mark_kmer);
    printf("Hashing subcontigs and finding unique k-mers\n");
    hash_and_insert(ht, subcontigs, hashtable_add_kmer);

    printf("A total of %ld different k-mers were found\n%ld k-mers were unique\n",ht->count, sum_unique_hahses(ht));

    // write tsv of unique hashes file
    FILE *kmercontent;
    char* subcontig_info;
    kmercontent = fopen(outdir, "w+");
    if(kmercontent == NULL){
        fprintf(stderr, "Error: failed to open the specified output directory, exiting\n");
        return EXIT_FAILURE;
    }
    fprintf(kmercontent,"SubcontigID\tStrainID\tContigID\tStart_Stop\tLength\tNunique\n");
    uint32_t i=1;
    while(ht->subcontig_names[i]!=NULL){
        fprintf(kmercontent,"%s\t", ht->subcontig_names[i]);
        subcontig_info = strtok(ht->subcontig_names[i], ";");
        for(int j=0; j<4; ++j){
            fprintf(kmercontent,"%s\t", subcontig_info);
            subcontig_info = strtok(NULL, ";");
        }
        fprintf(kmercontent,"%d\n", ht->subcontig_counts[i]);
        ++i;
    }

    printf("K-mers hashed and counted, the results can be found in the output directory under KmerContent.report\n");

    fclose(kmercontent);
    free(outdir);
    free(subcontigs);
    free(exc_subcontigs);
    hashtable_destroy(ht);
}