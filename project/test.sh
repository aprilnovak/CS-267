#!/bin/bash

# remove old timing file since we will append to a clean file
rm timing.txt
touch timing.txt

# run for a range of elements
for i in $(seq 1000 1000 10000)
do
  echo "Running for $i elements ..."
  ./a.out 4 $i 1 1 2
done
