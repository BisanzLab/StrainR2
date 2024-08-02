#include <dirent.h>
#include <getopt.h>
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
    "\t\t-h\t\t\t: display this message again\n"

typedef enum ht_element_status{
    EMPTY = 0,
    UNIQUE = 1,
    NON_UNIQUE = 2
} ht_element_status;

typedef struct ht_element{
    uint64_t key; // key is a hash
    ht_element_status status; // status of hash
    uint32_t subcontig_id; // id is array index for hash name
} ht_element;

typedef struct ht_element_small{
    uint32_t key; // key is a hash
    uint32_t value; // upper 4 bits are status, lower 28 are subcontig id
} ht_element_small;

typedef struct hashtable{
    ht_element* items;
    char** subcontig_names;
    uint64_t size;
    uint64_t entry_bitmask;
    uint64_t count;
    uint32_t* subcontig_counts;
    uint32_t num_subcontigs;
    uint32_t kmer_size;
    ht_element_small* items_small; // for use in memory-efficient option
    bool is_small;
} hashtable;


hashtable* hashtable_create(uint32_t kmer_size, bool is_small, uint32_t num_subconts);
void hashtable_destroy(hashtable* ht);ht_element* hashtable_insert(hashtable* ht, uint64_t key, ht_element_status status, uint32_t subcontig_id);
ht_element_small* hashtable_insert_small(hashtable* ht, uint32_t key, ht_element_status status, uint32_t subcontig_id);
void hashtable_resize(hashtable* ht);
void hashtable_resize_small(hashtable* ht);
uint64_t MurmurHash64A (const void* key, int len, uint64_t seed);
uint32_t MurmurHash3_x86_32(const void * key, int len, uint32_t seed);
void hash_and_insert_subcontig(hashtable* ht, char* seq, uint32_t subcontig_id, void (*kmer_func)(hashtable*, char*, uint32_t));
void hash_and_insert(hashtable* ht, char* dir_location, void (*kmer_func)(hashtable*, char*, uint32_t));
