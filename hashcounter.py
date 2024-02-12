#!/usr/bin/env python3
import sourmash
import screed
import os
import sys
import time
import csv


unique_hash_count = {}
#key = subcontig name
#value = count
unique_hashes = {}
#key = hash
#value = 
        #None if never seen kmer before
        #int for originating subcontig num if assumed unique
        #False is not unique
subcont_name_number = {}
#key = subcontig number (used to save memory)
#value = name


def hash_and_count():
        i=0
        print(f"Reading {format(len(excludedSubcontigs),',d')} excluded subcontigs (subcontigs below min subcontig size threshold) and logging hashes")
        for excludedSubcontig in excludedSubcontigs:
                i+=1
                #create hashing object
                mh = sourmash.MinHash(0, k, scaled=1)

                #load subcontig sequence from fasta and hash
                for seq in screed.open(outdir + "excludedSubcontigs/" + excludedSubcontig):
                        mh_list = mh.seq_to_hashes(seq.sequence, force=True)
                        name = seq.name
                        subcont_name_number[i]=name
                        unique_hash_count[name] = 0
                
                #consider all hashes from excluded subcontigs as non unique so that if they appear in an included subcontig they are not counted
                for mh in mh_list:
                        if mh not in unique_hashes:
                                unique_hashes[mh] = False
                del mh_list
                del mh
                del seq

        print(f"Found {format(len(unique_hashes),',d')} unique hashes in all excluded subcontigs -- marking them as non unique")
        print(f"Reading and hashing {format(len(subcontigs),',d')} non-excluded subcontigs - counting unique hashes")
        for subcontig in subcontigs:
                i+=1
                #create hashing object
                mh = sourmash.MinHash(0, k, scaled=1)

                #load subcontig sequence from fasta and hash
                for seq in screed.open(outdir + "Subcontigs/" + subcontig):
                        mh_list = mh.seq_to_hashes(seq.sequence, force=True)
                        name = seq.name
                        subcont_name_number[i]=name
                        unique_hash_count[name] = 0

                #go through all hashes and get counts of how many are unique per subcontig
                for mh in mh_list:

                        if mh not in unique_hashes:
                                unique_hash_count[name] += 1
                                unique_hashes[mh] = i
                        elif type(unique_hashes.get(mh)) is int:
                                unique_hash_count[subcont_name_number[unique_hashes[mh]]] -= 1
                                unique_hashes[mh] = False
                del mh_list
                del mh
                del seq

        print(f"Found total of {format(len(unique_hashes),',d')} unique hashes")
        print("writing unique hash counts to KmerContent.report")
        #write to tsv
        with open(outdir + "KmerContent.report", "w") as kmer_report:
                tsv_writer = csv.writer(kmer_report, delimiter="\t", lineterminator="\n")
                tsv_writer.writerow(["SubcontigID", "StrainID", "ContigID", "Start_Stop", "Length", "Nunique"])

                for SubcontigID in unique_hash_count:
                        subcontig_split = SubcontigID.split(";")
                        if(subcontig_split[1][0:9]!="EXCLUDED_"):
                                tsv_writer.writerow([SubcontigID, subcontig_split[0], subcontig_split[1], subcontig_split[2], subcontig_split[3], unique_hash_count[SubcontigID]])



if __name__ == "__main__":
        #proper arg parsing is done by PreProcessR script calling on HashCounter.py
        #outdir is first argument passed to HashCounter from PreprocessR, where subcontigs and .report file will be
        #kmersize is second
        if not len(sys.argv) == 3:
                print(f"Missing Inputs to Hash Counter! Got {sys.argv}\nFirst arg is path to db, Second arg is kmersize as int")
                raise SystemExit 
        
        outdir = str(sys.argv[1]) + "/"
        k = int(sys.argv[2])
        start_time = time.time()
        subcontigs = os.listdir(outdir + "Subcontigs")
        excludedSubcontigs = os.listdir(outdir + "excludedSubcontigs")

        print(f"Started counting unique hashes\nK-Size: {k}\nReference Directory: {outdir}")
        hash_and_count()
        print(f"Hashes counted in {int((time.time() - start_time)//60)} minutes and {int((time.time() - start_time)%60)} seconds")