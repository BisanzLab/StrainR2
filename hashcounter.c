#include "hashcounter.h"
#include "kseq.h"

/*
 * Implemented hashtable has open addressing with linear probe collision policy
 * A memory efficient hashtable entry is also available, with half the memory usage (collisions are more likely with this)
 * The hashtable resizes when the load factor exceeds 0.75 after entering the k-mers of a subcontig
 * Keys are k-mers hashed using the non-cryptographic MurMurHash
 * Values are the status of the k-mer (i.e. unique or not) and also the id of the subcontig from which it originates
 */

KSEQ_INIT(int, read);

hashtable* hashtable_create(uint32_t kmer_size, bool is_small, uint32_t num_subconts){
    hashtable* ht = (hashtable*) malloc(sizeof(hashtable));
    ht->subcontig_names = calloc(num_subconts,sizeof(char*));
    ht->subcontig_counts = calloc(num_subconts,sizeof(int));
    ht->num_subcontigs = num_subconts;
    ht->size = INITIAL_HT_SIZE;
    ht->count = 0;
    ht->entry_bitmask = INITIAL_HT_BITMASK;
    ht->kmer_size = kmer_size;
    ht->is_small = is_small;
    if(is_small){
        ht->items_small = (ht_element_small*) calloc(INITIAL_HT_SIZE, sizeof(ht_element_small));
    }else{
        ht->items = (ht_element*) calloc(INITIAL_HT_SIZE, sizeof(ht_element));
    }
    return ht;
}

void hashtable_destroy(hashtable* ht){
    int i=1;
    while(ht->subcontig_names[i]!=NULL){
        free(ht->subcontig_names[i]);
        ++i;
    }
    free(ht->subcontig_names);
    free(ht->subcontig_counts);
    if(ht->is_small){
        free(ht->items_small);
    }else{
        free(ht->items);
    }
    free(ht);
}

static inline ht_element_status ht_small_get_status(ht_element_small* element){
    return (element->value & 0xC0000000) >> 30;
}

static inline void ht_small_set_status(ht_element_small* element, ht_element_status status){
    element->value &= 0x3FFFFFFF;
    element->value |= (status << 30);
}

static inline uint32_t ht_small_get_id(ht_element_small* element){
    return element->value & 0x3FFFFFFF;
}

static inline void ht_small_set_id(ht_element_small* element, uint32_t id){
    element->value &= 0xC0000000;
    element->value |= id;
}

// insert into hashtable with linear probe collision policy
// return entry if found in ht 
ht_element* hashtable_insert(hashtable* ht, uint64_t key, ht_element_status status, uint32_t subcontig_id){
    uint64_t hash = key & ht->entry_bitmask;
    ht_element* current_item = &(ht->items[hash]);
    while(current_item->status != EMPTY){
        if(current_item->key == key) return current_item;
        ++current_item;
        // reset to beginning of hashtable if end is reached
        if(current_item - ht->items == ht->size) current_item = ht->items;
    }
    current_item->key = key;
    current_item->status = status;
    current_item->subcontig_id = subcontig_id;
    return NULL;
}

// insert function with same behaviour, but for memory-efficient version
ht_element_small* hashtable_insert_small(hashtable* ht, uint32_t key, ht_element_status status, uint32_t subcontig_id){
    uint32_t hash = key & ht->entry_bitmask;
    ht_element_small* current_item = &(ht->items_small[hash]);
    while(ht_small_get_status(current_item) != EMPTY){
        if(current_item->key == key) return current_item;
        ++current_item;
        // reset to beginning of hashtable if end is reached
        if(current_item - ht->items_small == ht->size) current_item = ht->items_small;
    }
    current_item->key = key;
    ht_small_set_status(current_item, status);
    ht_small_set_id(current_item, subcontig_id);
    return NULL;
}

// double ht size and re-enter all elements from left to right
void hashtable_resize(hashtable* ht){
    if(ht->is_small) return hashtable_resize_small(ht);
    ht->size *= 2;
    printf("Hashtable is resizing, new size will use ~ %ld GiB of memory\n", ht->size/INITIAL_HT_SIZE/2);
    uint64_t changed_bit = ht->entry_bitmask;
    ht->entry_bitmask = (ht->entry_bitmask << 1) | 0x1;
    changed_bit ^= ht->entry_bitmask;
    ht->items = realloc(ht->items, ht->size * sizeof(ht_element));
    ht_element* current_entry =  ht->items-1;
    while(current_entry != &ht->items[ht->size/2]){
        ++current_entry;
        if(current_entry->status == EMPTY) continue;
        if(hashtable_insert(ht, current_entry->key, current_entry->status, current_entry->subcontig_id)!=NULL) continue;
        current_entry->key = 0;
        current_entry->status = EMPTY;
        current_entry->subcontig_id = 0;
    }
}

void hashtable_resize_small(hashtable* ht){
    /*
    add warnings and errors for big sizes
    */
    ht->size *= 2;
    printf("Hashtable is resizing, new size will use ~ %.1f GiB of memory\n", (float)ht->size/INITIAL_HT_SIZE/4);
    if(ht->size == 536870912){ // 2^29
        printf("Warning: Due to the large input size and use of the memory-efficient mode, the output is losing some accuracy.\n");
    }
    if(ht->size == 1073741824){ // 2^30
        fprintf(stderr,"Error: memory-efficient mode has lost too much accuracy to continue, please try again without it enabled.\n");
        exit(EXIT_FAILURE);
    }
    uint64_t changed_bit = ht->entry_bitmask;
    ht->entry_bitmask = (ht->entry_bitmask << 1) | 0x1;
    changed_bit ^= ht->entry_bitmask;
    ht->items_small = realloc(ht->items_small, ht->size * sizeof(ht_element_small));
    ht_element_small* current_entry =  ht->items_small-1;
    while(current_entry != &ht->items_small[ht->size/2]){
        ++current_entry;
        if(ht_small_get_status(current_entry) == EMPTY) continue;
        if(hashtable_insert_small(ht, current_entry->key, ht_small_get_status(current_entry), ht_small_get_id(current_entry))!=NULL) continue;
        current_entry->key = 0;
        ht_small_set_status(current_entry, EMPTY);
        ht_small_set_id(current_entry, 0);
    }
}

// return the sum of all unique hashes
uint64_t sum_unique_hahses(hashtable* ht){
    uint64_t sum = 0;
    for(int i=0; i<ht->num_subcontigs; ++i){
        sum+=ht->subcontig_counts[i];
    }
    return sum;
}

// add one k-mer (to be hashed) to the hashtable
static inline void hashtable_add_kmer(hashtable* ht, char* seq, uint32_t subcontig_id){
    uint64_t hash = MurmurHash64A(seq, ht->kmer_size, (uint64_t)07062024);
    ht_element* hashtable_item = hashtable_insert(ht, hash, UNIQUE, subcontig_id);
    if(hashtable_item == NULL){
        ++ht->subcontig_counts[subcontig_id];
        ++ht->count;
    } else if(hashtable_item->status == UNIQUE){
        hashtable_item->status = NON_UNIQUE;
        --ht->subcontig_counts[hashtable_item->subcontig_id];
    }
}

// same functionality but for memory-efficient mode
static inline void hashtable_small_add_kmer(hashtable* ht, char* seq, uint32_t subcontig_id){
    uint32_t hash = MurmurHash3_x86_32(seq, ht->kmer_size, (uint32_t)07062024);
    ht_element_small* hashtable_item = hashtable_insert_small(ht, hash, UNIQUE, subcontig_id);
    if(hashtable_item == NULL){
        ++ht->subcontig_counts[subcontig_id];
        ++ht->count;
    } else if(ht_small_get_status(hashtable_item) == UNIQUE){
        ht_small_set_status(hashtable_item, NON_UNIQUE);
        --ht->subcontig_counts[ht_small_get_id(hashtable_item)];
    }
}

// function for marking a k-mer as non-unique
static inline void hashtable_mark_kmer(hashtable* ht, char* seq, uint32_t subcont_id){
    uint64_t hash = MurmurHash64A(seq, ht->kmer_size, (uint64_t)07062024);
    if(hashtable_insert(ht, hash, NON_UNIQUE, subcont_id)==NULL) ++ht->count;
}

static inline void hashtable_small_mark_kmer(hashtable* ht, char* seq, uint32_t subcont_id){
    uint32_t hash = MurmurHash3_x86_32(seq, ht->kmer_size, (uint32_t)07062024);
    if(hashtable_insert_small(ht, hash, NON_UNIQUE, subcont_id)==NULL) ++ht->count;
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
        case 7: h ^= (uint64_t)(data2[6]) << 48;
        case 6: h ^= (uint64_t)(data2[5]) << 40;
        case 5: h ^= (uint64_t)(data2[4]) << 32;
        case 4: h ^= (uint64_t)(data2[3]) << 24;
        case 3: h ^= (uint64_t)(data2[2]) << 16;
        case 2: h ^= (uint64_t)(data2[1]) << 8;
        case 1: h ^= (uint64_t)(data2[0]);
            h *= m;
    };
    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

/* following function adapted from Austin Appleby */
uint32_t MurmurHash3_x86_32(const void* key, int len, uint32_t seed){
    const uint8_t* data = (const uint8_t*)key;
    const int nblocks = len / 4;
    uint32_t h1 = seed;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    const uint32_t* blocks = (const uint32_t*)(data + nblocks*4);
    for(int i = -nblocks; i; ++i){
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> 17);
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> 19); 
        h1 = h1*5+0xe6546b64;
    }

    const uint8_t* tail = (const uint8_t*)(data + nblocks*4);
    uint32_t k1 = 0;
    switch(len & 3){
    case 3: k1 ^= tail[2] << 16;
    case 2: k1 ^= tail[1] << 8;
    case 1: k1 ^= tail[0];
            k1 *= c1; k1 = (k1 << 15) | (k1 >> 17); k1 *= c2; h1 ^= k1;
    };

    h1 ^= len;
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
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

// add k-mers to the hashtable for an entire subcontig
void hash_and_insert_subcontig(hashtable* ht, char* seq, uint32_t subcontig_id, void (*kmer_func)(hashtable*, char*, uint32_t)){
    uint32_t i = 0;
    char* rc = reverse_complement(seq);
    uint32_t n;
    uint32_t seq_len = strlen(seq);
    while(n=check_n(&seq[i], ht->kmer_size)){
        i+=n;
    }
    while(ht->kmer_size + i <= seq_len){
        if(seq[i+ht->kmer_size-1]=='N'){
            i+=ht->kmer_size;
            if(i+ht->kmer_size > seq_len) break;
            while(n=check_n(&seq[i], ht->kmer_size)) i+=n;
            continue;
        }
    
        if(strncmp(&seq[i], &rc[seq_len-ht->kmer_size-i], ht->kmer_size) < 0){
            kmer_func(ht, &seq[i], subcontig_id);
        }else{
            kmer_func(ht, &rc[seq_len-ht->kmer_size-i], subcontig_id);
        }
        
        ++i;
    }
    free(rc);
    // resize hashtable if load factor is >0.75 after subcontig addition
    if((float) ht->count / ht->size > 0.75) hashtable_resize(ht);
}

// add k-mers to the hashtable for all subcontigs in a directory
void hash_and_insert(hashtable* ht, char* dir_location, void (*kmer_func)(hashtable*, char*, uint32_t)){
    kseq_t* seq;
    struct dirent *de;
    DIR *dr = opendir(dir_location);
    if(dr == NULL) {
        fprintf(stderr, "Could not open excluded subcontigs directory\n\n");
        exit(EXIT_FAILURE);
    }
    char* subcont_location;
    char* subcont_name;
    uint32_t subcont_id;
    while (((de = readdir(dr)) != NULL)) {
        if(!(strlen(de->d_name) >= 10 && strcmp(&de->d_name[strlen(de->d_name) - 10], ".subcontig") == 0)) continue;
        uint32_t loc_size = strlen(dir_location)+strlen(de->d_name)+1;
        subcont_location = calloc(loc_size, sizeof(char));
        strcpy(subcont_location, dir_location);
        subcont_location = strcat(subcont_location, de->d_name);
        FILE* fp = fopen(subcont_location,"r");
        if(fp == NULL){
            fprintf(stderr, "Error opening %s\n", de->d_name);
            exit(EXIT_FAILURE);
        }
        seq = kseq_init(fileno(fp));
        kseq_read(seq);
        //printf("name: %s\n", seq->name.s);
        //printf("comment: %s\n", seq->comment.s);
        memcpy(subcont_location, de->d_name, strlen(de->d_name)+1);
        strtok(subcont_location, ".");
        strtok(subcont_location, "_");
        subcont_id = atoi(strtok(NULL,"_"));
        if(seq->comment.s != NULL){
            subcont_name = calloc(strlen(seq->name.s)+strlen(seq->comment.s)+1, sizeof(char));
            memcpy(subcont_name, seq->name.s, strlen(seq->name.s));
            strcat(subcont_name, seq->comment.s);
        } else {
            subcont_name = calloc(strlen(seq->name.s)+1, sizeof(char));
            memcpy(subcont_name, seq->name.s, strlen(seq->name.s));
        }
        ht->subcontig_names[subcont_id] = subcont_name;
        hash_and_insert_subcontig(ht, seq->seq.s, subcont_id, kmer_func);
        free(subcont_location);
        fclose(fp);
        kseq_destroy(seq);
    }
    closedir(dr);
}

int main(int argc, char **argv){
    int opt;
    char *subcontigs = NULL;
    char *exc_subcontigs = NULL;
    char *outdir = NULL;
    uint32_t kmer_size = 0;
    bool is_mem_efficient = false;
    uint32_t num_subcontigs = 0;

    // parse options
    while ((opt = getopt(argc, argv, "s:e:k:n:o:mh")) != -1) {
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
            case 'o': {
                outdir = calloc(strlen(optarg) + strlen("/kmercontent.report") + 1, sizeof(char));
                strcpy(outdir, optarg);
                strcat(outdir, "/kmercontent.report");
            } break;
            case 'm': {
                is_mem_efficient = true;
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

    if(is_mem_efficient){
        printf("Memory-efficient mode has been enabled. Note that this comes with reduced accuracy when there are larger input sizes.\n");
    }

    printf("Hashing and counting k-mers\n");
    hashtable* ht = hashtable_create(kmer_size, is_mem_efficient, num_subcontigs);

    // main pipeline
    if(is_mem_efficient){
        printf("Hashing excluded subcontigs and marking them as non-unique\n");
        hash_and_insert(ht, exc_subcontigs, hashtable_small_mark_kmer);
        printf("Hashing subcontigs and finding unique k-mers\n");
        hash_and_insert(ht, subcontigs, hashtable_small_add_kmer);
    }else{
        printf("Hashing excluded subcontigs and marking them as non-unique\n");
        hash_and_insert(ht, exc_subcontigs, hashtable_mark_kmer);
        printf("Hashing subcontigs and finding unique k-mers\n");
        hash_and_insert(ht, subcontigs, hashtable_add_kmer);
    }

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

    printf("K-mers hashed and counted, the results can be found in the output directory under kmercontent.report\n");

    fclose(kmercontent);
    free(outdir);
    free(subcontigs);
    free(exc_subcontigs);
    hashtable_destroy(ht);
}