#include <dirent.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define OVERLAP_LENGTH 500
#define USAGE                                                                                                                                        \
    "USAGE: subcontig -i path/to/in [OPTIONS]\n"                                                                                                     \
    "subcontig splits input genomes into smaller parts and save to .subcontig files\n"                                                               \
    "\tRequired Arguments:\n"                                                                                                                        \
    "\t\t-i path/to/genomes\t: path to the directory for all community genomes\n"                                                                    \
    "\tOptional Arguments:\n"                                                                                                                        \
    "\t\t-o path/to/out\t: path to your output directory [Default = current directory]\n"                                                            \
    "\t\t-s number\t: maximum subcontig size (overrides default use of smallest N50) [Default = calculated N50]\n"                                   \
    "\t\t-e number\t: exclude subcontig size (minimum subcontig size) [Default = 10000]\n"                                                           \
    "\t\t-h\t\t: display this message again\n"

// for testing purposes this code was compiled with:
// gcc -g -fsanitize=address -std=gnu99 -Wall -Wextra -Werror -Wno-unused-function -Wno-unused-parameter -O0 -o subcontig subcontig.c

// the provided binary was compiled using:
// gcc -o subcontig subcontig.c

int subcontigCount = 0;

// save a sequence and appropriate header information to outdir (or excludedSubcontigs if it is less than minSubcontigSize)
void saveSubcontig(char *outdir, char *subcontigName, char *strainID, char *subcontigSeq, int start, int length, char *overlap);
// subcontig a genome and write the sequences to outdir
void writeSubcontigs(char *outdir, char *excludeDir, char *genomeLocation, char *strainID, int *contigLengths, int maxSubcontigSize,
                     int minSubcontigSize);
// returns array of contig lengths for a given genome and passes array size to contigLengthsSize
int *getContigLengths(char *genomeLocation, int minSubcontigSize, int *contigLengthsSize);
// compare function for qsort
int compare(const void *a, const void *b);

int main(int argc, char **argv) {

    // define options and their defaults
    int opt;
    char *indir = NULL;
    char *outdir = malloc(2 * sizeof(char));
    strcpy(outdir, ".");
    int minSubcontigSize = 10000;
    int maxSubcontigSize = 0;

    // parse options
    while ((opt = getopt(argc, argv, "i:o:e:s:h")) != -1) {
        switch (opt) {
            case 'i': {
                indir = calloc(strlen(optarg) + 1, sizeof(char));
                strcpy(indir, optarg);
            } break;
            case 'o': {
                outdir = realloc(outdir, (strlen(optarg) + 1) * sizeof(char));
                strcpy(outdir, optarg);
            } break;
            case 'e': {
                minSubcontigSize = atoi(optarg);
            } break;
            case 's': {
                maxSubcontigSize = atoi(optarg);
            } break;
            case 'h': {
                printf(USAGE);
                return EXIT_SUCCESS;
            }
            default: {
                printf(USAGE);
                free(outdir);
                return EXIT_FAILURE;
            }
        }
    }

    // check validity of inputs

    // check indir is given
    if (indir == NULL) {
        printf(USAGE);
        free(outdir);
        return EXIT_FAILURE;
    }

    // make output locations
    char *excludeDir;
    char *temp;
    excludeDir = calloc((strlen(outdir) + strlen("/excludedSubcontigs/") + 1), sizeof(char));
    sprintf(excludeDir, "%s/excludedSubcontigs/", outdir);

    temp = malloc((strlen(outdir)+1) * sizeof(char));
    strcpy(temp,outdir);
    outdir = realloc(outdir, (strlen(outdir) + strlen("/Subcontigs/") + 1) * sizeof(char));
    sprintf(outdir, "%s/Subcontigs/", temp);

    temp = realloc(temp, (strlen(indir)+1) * sizeof(char));
    strcpy(temp,indir);
    indir = realloc(indir, (strlen(indir) + 2) * sizeof(char));
    sprintf(indir, "%s/", temp);

    free(temp);


    // check outdirs don't already exist
    DIR *outd;
    DIR *excld;
    if ((outd = opendir(outdir)) != NULL || (excld = opendir(excludeDir)) != NULL) {
        fprintf(stderr, "Error: Subcontig files already exist in out-directory\n");
        if (outd != NULL) {
            closedir(outd);
        } else if (excld != NULL) {
            closedir(excld);
        }
        free(indir);
        free(outdir);
        free(excludeDir);
        return EXIT_FAILURE;
    }

    mkdir(outdir, 0777);
    mkdir(excludeDir, 0777);

    // check indir exists
    struct dirent *de;
    DIR *dr = opendir(indir);
    if (dr == NULL) {
        fprintf(stderr, "Could not open input directory\n\n");
        return EXIT_FAILURE;
    }

    if(minSubcontigSize >= 250000){
        fprintf(stderr, "Error: Minimum subcontig size is too large, please set it to be less than 250,000");
        return EXIT_FAILURE;
    }
    if(minSubcontigSize <= 5000){
        fprintf(stderr, "Warning: It is strongly recommended not to set Minimum subcontig size to be smaller than 5000 to ensure multi-copy elements are not present");
    }

    if(maxSubcontigSize==0){
        // calculate lowest N50
        int smallestN50 = 0;
        char* smallestN50_genome = malloc(0);
        while (((de = readdir(dr)) != NULL)) {
            if ((strlen(de->d_name) >= 4 && strcmp(&de->d_name[strlen(de->d_name) - 4], ".fna") == 0) ||
                (strlen(de->d_name) >= 6 && strcmp(&de->d_name[strlen(de->d_name) - 6], ".fasta") == 0)) {

                char *genomeLocation;
                genomeLocation = calloc(strlen(indir) + strlen(de->d_name) + 1, sizeof(char));
                sprintf(genomeLocation, "%s%s", indir, de->d_name);

                int numContigs = 0;
                int *contigLengths = getContigLengths(genomeLocation, minSubcontigSize, &numContigs);
                // find N50
                qsort(contigLengths, numContigs, sizeof(int), compare);
                int sum = 0;
                for (int i = 0; i < numContigs; ++i) {
                    if(contigLengths[i] > minSubcontigSize){
                        sum += contigLengths[i];
                    }
                }
                int i = 0;
                int contigSum = 0;
                while (contigSum < sum / 2) {
                    if(contigLengths[i] > minSubcontigSize){
                        contigSum += contigLengths[i];
                    }
                    ++i;
                }
                int N50 = i!=0 ? contigLengths[i-1] : contigLengths[0];
                // see if that N50 is the smallest one
                if(N50 < smallestN50 || smallestN50 == 0){
                    smallestN50 = N50;
                    smallestN50_genome = realloc(smallestN50_genome, sizeof(char)*strlen(genomeLocation)+1);
                    strcpy(smallestN50_genome, genomeLocation);
                }

                free(contigLengths);
                free(genomeLocation);
            }
        }
        maxSubcontigSize = smallestN50;
        printf("Smallest N50 is %d, which belongs to %s\n", maxSubcontigSize, smallestN50_genome);
        free(smallestN50_genome);
    }

    if(maxSubcontigSize > 250000){
        fprintf(stdout, "Warning: Subcontig size can not be larger than 250Kb -- setting subcontig size to 250,000");
        maxSubcontigSize = 250000;
    }
    if(maxSubcontigSize < minSubcontigSize + 1000){
        fprintf(stdout, "Warning: Max subcontig size is too small -- setting it to the exclude size + 1000 bases (%d bases)", minSubcontigSize + 1000);
        maxSubcontigSize = minSubcontigSize + 1000;
    }

    closedir(dr);
    dr = opendir(indir);
    if (dr == NULL) {
        fprintf(stderr, "Could not open input directory\n\n");
        return EXIT_FAILURE;
    }

    // write subcontigs
    while (((de = readdir(dr)) != NULL)) {
        if ((strlen(de->d_name) >= 4 && strcmp(&de->d_name[strlen(de->d_name) - 4], ".fna") == 0) ||
            (strlen(de->d_name) >= 6 && strcmp(&de->d_name[strlen(de->d_name) - 6], ".fasta") == 0)) {

            char *genomeLocation;
            genomeLocation = calloc(strlen(indir) + strlen(de->d_name) + 1, sizeof(char));
            sprintf(genomeLocation, "%s%s", indir, de->d_name);

            int *contigLengths = getContigLengths(genomeLocation, minSubcontigSize, NULL);
            writeSubcontigs(outdir, excludeDir, genomeLocation, strtok(de->d_name, "."), contigLengths, maxSubcontigSize, minSubcontigSize);
            free(contigLengths);
            free(genomeLocation);
        }
    }
    closedir(dr);
    free(indir);
    free(outdir);
    free(excludeDir);

    return EXIT_SUCCESS;
}

// split contigs into subcontigs and write the sequences to outdir
void writeSubcontigs(char *outdir, char *excludeDir, char *genomeLocation, char *strainID, int *contigLengths, int maxSubcontigSize,
                     int minSubcontigSize) {
    FILE *genome = fopen(genomeLocation, "r");
    char *line = NULL;
    size_t maxLen = 0;
    ssize_t lineLen = -1;
    char *subcontigName;
    char *subcontigSeq;
    char *overlapBuff;
    int seqIndex = 0;
    int contigIndex = 0;
    int subcontigLengths = 0;
    subcontigLengths = contigLengths[contigIndex] / (contigLengths[contigIndex] / (maxSubcontigSize+1) + 1);
    int start = 1;

    if (genome == NULL) {
        fprintf(stderr, "Error opening %s\n\n", genomeLocation);
        exit(EXIT_FAILURE);
    }

    lineLen = getline(&line, &maxLen, genome);
    subcontigName = calloc(lineLen - 1, sizeof(char));
    strncpy(subcontigName, &line[1], lineLen - 2);
    subcontigSeq = calloc(subcontigLengths + 1, sizeof(char));
    overlapBuff = calloc(OVERLAP_LENGTH + 1, sizeof(char));

    // read genome files line by line 
    while ((lineLen = getline(&line, &maxLen, genome)) != -1) {

        if (line[0] == '>') {
            free(overlapBuff);
            overlapBuff = NULL;
            if (seqIndex >= minSubcontigSize) {
                saveSubcontig(outdir, subcontigName, strainID, subcontigSeq, start, seqIndex, overlapBuff);
            } else if (seqIndex >= OVERLAP_LENGTH) {
                char *excludedSubcontigName = calloc(strlen(subcontigName) + 10, sizeof(char));
                sprintf(excludedSubcontigName, "EXCLUDED_%s", subcontigName);
                saveSubcontig(excludeDir, excludedSubcontigName, strainID, subcontigSeq, start, seqIndex, overlapBuff);
                free(excludedSubcontigName);
            }
            overlapBuff = calloc(OVERLAP_LENGTH + 1, sizeof(char));
            start += seqIndex;
            ++contigIndex;
            seqIndex = 0;
            free(subcontigSeq);
            free(subcontigName);
            subcontigLengths = contigLengths[contigIndex] / (contigLengths[contigIndex] / (maxSubcontigSize+1) + 1);
            subcontigSeq = calloc(subcontigLengths + 1, sizeof(char));
            subcontigName = calloc(lineLen - 1, sizeof(char));
            strncpy(subcontigName, &line[1], lineLen - 2);

        } else if (seqIndex + lineLen - 1 <= subcontigLengths) {
            strncpy(&subcontigSeq[seqIndex], line, lineLen - 1);
            seqIndex += lineLen - 1;
        } else {
            strncpy(&subcontigSeq[seqIndex], line, subcontigLengths - seqIndex);
            saveSubcontig(outdir, subcontigName, strainID, subcontigSeq, start, subcontigLengths, overlapBuff);
            strcpy(overlapBuff, &subcontigSeq[strlen(subcontigSeq) - OVERLAP_LENGTH]);
            free(subcontigSeq);

            // handle the case where a line is larger than the subcontig size
            int i = 1;
            while(lineLen - (subcontigLengths*i - seqIndex) > subcontigLengths){
                subcontigSeq = calloc(subcontigLengths + 1, sizeof(char));
                strncpy(subcontigSeq, &line[subcontigLengths*i + seqIndex], subcontigLengths);
                start += subcontigLengths;
                saveSubcontig(outdir, subcontigName, strainID, subcontigSeq, start, subcontigLengths, overlapBuff);
                strcpy(overlapBuff, &subcontigSeq[strlen(subcontigSeq) - OVERLAP_LENGTH]);
                free(subcontigSeq);
                ++i;
            }

            subcontigSeq = calloc(subcontigLengths + 1, sizeof(char));
            char* resized_line = calloc(strlen(line),sizeof(char));
            strncpy(resized_line, line, strlen(line)-1);
            strcpy(subcontigSeq, &resized_line[subcontigLengths*i - seqIndex]);
            start += subcontigLengths;
            seqIndex = strlen(subcontigSeq);
        }
    }

    free(subcontigSeq);
    free(overlapBuff);
    free(line);
    free(subcontigName);
    fclose(genome);
}

// save a sequence and appropriate header information to outdir (or excludedSubcontigs if it is less than minSubcontigSize)
// (called from by writeSubcontigs)
void saveSubcontig(char *outdir, char *subcontigName, char *strainID, char *subcontigSeq, int start, int length, char *overlap) {
    FILE *fptr = NULL;
    ++subcontigCount;
    char *seq = NULL;
    int overlapLen = 0;

    if (overlap == NULL) {
        seq = calloc(length + 1, sizeof(char));
        strcpy(seq, subcontigSeq);
    } else {
        seq = calloc(length + strlen(overlap) + 1, sizeof(char));
        strcpy(seq, overlap);
        strcpy(&seq[strlen(overlap)], subcontigSeq);
        overlapLen = strlen(overlap);
    }

    size_t needed = snprintf(NULL, 0, ">%s;%s;%d_%d;%d", strainID, subcontigName, start - overlapLen, start + length - 1, length + overlapLen) + 1;
    char *savedSubcontigName = calloc(needed, sizeof(char));
    sprintf(savedSubcontigName, ">%s;%s;%d_%d;%d", strainID, subcontigName, start - overlapLen, start + length - 1, length + overlapLen);

    needed = snprintf(NULL, 0, "%ssubcontig_%d.subcontig", outdir, subcontigCount) + 1;
    char *subcontigLocation = calloc(needed, sizeof(char));
    sprintf(subcontigLocation, "%ssubcontig_%d.subcontig", outdir, subcontigCount);

    fptr = fopen(subcontigLocation, "w");
    if (fptr == NULL) {
        fprintf(stderr, "Error writing %s\n\n", subcontigLocation);
        exit(EXIT_FAILURE);
    }

    fprintf(fptr, "%s\n", savedSubcontigName);
    int subcontigIndex = 0;
    char *newLine;
    newLine = calloc(81, sizeof(char));

    // write subcontig files line by line
    while (subcontigIndex < length + overlapLen) {
        strncpy(newLine, &seq[subcontigIndex], 80);
        subcontigIndex += 80;
        fprintf(fptr, "%s\n", newLine);
    }

    free(newLine);
    free(subcontigLocation);
    free(savedSubcontigName);
    free(seq);
    fclose(fptr);
}

// returns array of contig lengths for a given genome and passes array size to contigLengthsSize
int *getContigLengths(char *genomeLocation, int minSubcontigSize, int *contigLengthsSize) {
    FILE *genome = fopen(genomeLocation, "r");
    ssize_t lineLen = -1;
    size_t maxLen = 0;
    int contigIndex = 0;
    char *line = NULL;
    int *contigLengths = calloc(100, sizeof(int));
    int maxContigs = 100;
    int contigLength = 0;

    if (genome == NULL) {
        fprintf(stderr, "Error opening %s\n\n", genomeLocation);
        exit(EXIT_FAILURE);
    }

    lineLen = getline(&line, &maxLen, genome);

    while ((lineLen = getline(&line, &maxLen, genome)) != -1) {
        if (line[0] == '>') {
            if (maxContigs == contigIndex) {
                contigLengths = realloc(contigLengths, (maxContigs + 100) * sizeof(int));
                maxContigs += 100;
            }
            contigLengths[contigIndex] = contigLength;
            ++contigIndex;
            contigLength = 0;
        } else {
            contigLength += lineLen - 1;
        }
    }
    contigLengths[contigIndex] = contigLength;
    contigLengths = realloc(contigLengths, (contigIndex + 1) * sizeof(int));

    free(line);
    fclose(genome);
    if (contigLengthsSize != NULL) {
        *contigLengthsSize = contigIndex + 1;
    }
    return contigLengths;
}

// compare function used in sorting subcontig sizes and finding N50
int compare(const void *a, const void *b) {
    int *x = (int *)a;
    int *y = (int *)b;
    return *x - *y;
}