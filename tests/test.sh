#! /bin/bash

set -e

printf "Starting testing\n"
for test in ../tests/genomes/*/; do
  test_name=$(sed -E 's|.*/(.+)/$|\1|' <(echo $test))
  printf "\nTesting %s case\nSubcontig:\n" $test_name
  # subcontig testing
  ../src/subcontig -i $test -o ../tests
  diff -q ../tests/Subcontigs ../tests/expected_output/Subcontigs_"$test_name"
  diff -q ../tests/excludedSubcontigs ../tests/expected_output/excludedSubcontigs_"$test_name"
  rm -r ../tests/excludedSubcontigs ../tests/Subcontigs
  # hashcounter testing
  printf "Hashcounter:\n"
  ../src/hashcounter -s ../tests/expected_output/Subcontigs_$test_name \
    -e ../tests/expected_output/excludedSubcontigs_$test_name -o ../tests \
    -k 301 -n $(printf "%s+%s\n" $(ls -l ../tests/expected_output/Subcontigs_"$test_name" | wc -l) $(ls -l ../tests/expected_output/excludedSubcontigs_"$test_name" | wc -l) | bc)
  diff -q ../tests/KmerContent.report ../tests/expected_output/KmerContent_"$test_name".report
  rm  ../tests/KmerContent.report
done
