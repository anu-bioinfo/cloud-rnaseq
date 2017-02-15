#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: $0 <GDS_ID> <NUM_PARTITIONS>"
  exit 1
fi


DOC_ID=$1
NUM_PA=$2

p=0
while [ $p -lt $NUM_PA ]
do
  echo aegea batch submit --execute sra_download.py --storage /mnt=500 --ecr-image sra_download --environment GDS_ID=$DOC_ID NUM_PARTITIONS=$NUM_PA PARTITION_ID=$p RDB_S3_PATH=s3://czi-hca/config/dump.rdb.gz --memory 12000
  echo sleep 60
  p=$(( $p + 1))
done
