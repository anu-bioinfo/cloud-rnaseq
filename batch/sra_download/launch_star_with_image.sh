#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: $0 <GDS_IDS> <NUM_PARTITIONS>"
  exit 1
fi


DOC_IDS=$1
NUM_PA=$2

p=0
while [ $p -lt $NUM_PA ]
do
  echo aegea batch submit --execute run_star_and_htseq.py --storage /mnt=500 --ecr-image sra_download --environment GDS_IDS=$DOC_IDS NUM_PARTITIONS=$NUM_PA PARTITION_ID=$p --memory 64000
  echo sleep 100
  p=$(( $p + 1))
done
