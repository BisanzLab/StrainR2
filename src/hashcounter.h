#include <dirent.h>
#include <pthread.h>
#include <getopt.h>
#include <semaphore.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#define INITIAL_HT_SIZE 33554432 // 2^25 entries, hashtable will initially use 0.5 GiB in memory
#define INITIAL_HT_BITMASK 0x1FFFFFF // 25 1s
#define OVERLAP_LENGTH 500
#define USAGE                                                                                                                                        \
    "USAGE: hashcounter -s path/to/subconts -e path/to/exc_subconts -k kmer_size -o path/to/outdir\n"                                                \
    "hashcounter creates a log of how many kmers are unique in each subcontig, with excluded subcontig kmers considered non-unique\n"                \
    "\tRequired Arguments:\n"                                                                                                                        \
    "\t\t-s path/to/subconts\t: path to the directory for all subcontigs for which a report will be created\n"                                       \
    "\t\t-e path/to/exc_subconts\t: path to subcontigs whose kmers shall be considered non-unique\n"                                                 \
    "\t\t-k number\t\t: kmer sizes to use\n"                                                                                                         \
    "\t\t-n number\t\t: number of subcontigs (excluded or not) that will be input\n"                                                                 \
    "\t\t-o path/to/outdir\t: Directory to write output file to\n"                                                                                   \
    "\tOptional Arguments:\n"                                                                                                                        \
    "\t\t-t number\t: Number of threads to use [Default = 8]\n"                                                                                      \
    "\t\t-h\t\t\t: display this message again\n"

typedef enum ht_el_status{
    EMPTY = 0,
    UNIQUE = 1,
    NON_UNIQUE = 2
} ht_element_status;

typedef struct ht_element_t{
    uint64_t key; // key is a hash
    ht_element_status status; // status of hash
    uint32_t subcontig_id; // id is array index for hash name
} ht_element_t;

typedef struct hashtable_t{
    ht_element_t* items;
    char** subcontig_names;
    uint64_t size;
    uint64_t entry_bitmask;
    _Atomic uint64_t count;
    _Atomic int32_t* subcontig_counts;
    uint32_t num_subcontigs;
    uint32_t curr_subcontig; // non-atomic; only access in parent thread
    uint32_t kmer_size;
    uint32_t num_threads;
    sem_t can_edit_ht;
    pthread_mutex_t resize_lock;
} hashtable_t;

typedef void (*kmer_func_t)(hashtable_t*, char*, uint32_t, int32_t*, uint64_t*);

typedef struct hash_job_t{
    hashtable_t* ht;
    char* dir_location;
    char* subcontig_loc;
    uint64_t subcontig_id;
    kmer_func_t kmer_func;
    int32_t* local_uniq_counts;
} hash_job_t;

typedef struct job_queue_t{
    hash_job_t* job_head;
    pthread_mutex_t job_queue_lock;
} job_queue_t;

hashtable_t* hashtable_create(uint32_t kmer_size, uint32_t num_subconts, uint32_t num_threads);
void hashtable_destroy(hashtable_t* ht);
ht_element_t* hashtable_insert(hashtable_t* ht, uint64_t key, ht_element_status status, uint32_t subcontig_id);
void hashtable_resize(hashtable_t* ht);
uint64_t MurmurHash64A (const void* key, int len, uint64_t seed);
void* start_hash_thread(void* input);
void hash_and_insert_subcontig(hash_job_t* hash_job);
void hash_and_insert(hashtable_t* ht, char* dir_location, kmer_func_t kmer_func);
