#! /bin/bash

set -e

printf "Starting testing\n"
for test in ../tests/genomes/*/; do
  test_name=$(sed -E 's|.*/(.+)/$|\1|' <(echo $test))
  printf "\nTesting %s case\nSubcontig:\n" $test_name
  # subcontig testing
  ../src/subcontig -i $test -o ../tests
  if [ ! -z "$( ls -A ../tests/excludedSubcontigs )" ]; then 
    diff ../tests/excludedSubcontigs ../tests/expected_output/excludedSubcontigs_"$test_name"
  fi
  diff ../tests/Subcontigs ../tests/expected_output/Subcontigs_"$test_name"
  # hashcounter testing
  printf "Hashcounter:\n"
  ../src/hashcounter -s ../tests/Subcontigs -e ../tests/excludedSubcontigs -o ../tests \
    -k 301 -n $(printf "%s+%s\n" $(ls -l ../tests/Subcontigs | wc -l) $(ls -l ../tests/excludedSubcontigs | wc -l) | bc)
  diff <(sort ../tests/KmerContent.report) <(sort ../tests/expected_output/KmerContent_"$test_name".report)
  rm -r ../tests/excludedSubcontigs ../tests/Subcontigs
  rm  ../tests/KmerContent.report
done

printf "Testing successful\n"