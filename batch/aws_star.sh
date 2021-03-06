#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Usage: $0 <S3_BUCKET> <EXP_IDS> <NUM_PARTITIONS> <TAXON> <S3_SCRIPT_LOC>"
  echo "   <S3_BUCKET>: S3 Path. Example: s3://czi-hca/data"
  echo "   <EXP_IDS>: experiment directories under <S3_BUCKET> seperated by ','. Example: 200067835,200057872"
  echo "   <NUM_PARTITIONS>: Number of partitions you want to run. Example: 10 means you will run in 10 different dockercontainers"
  echo "   <TAXON>: mus or homo. mus will use mm10 as reference genome, homo will use hg38 as reference genome"
  echo "   <S3_SCRIPT_LOC>: script location. i.e. s3://czi-hca/config/run_star_and_htseq.py"
  exit 1
fi

S3_BUCKET=$1
DOC_IDS=$2
NUM_PA=$3
TAXON=$4
S3_SCRIPT_LOC=$5
SCRIPT_NAME=`basename $S3_SCRIPT_LOC`


p=0
while [ $p -lt $NUM_PA ]
do
  COMMAND="aws s3 cp $S3_SCRIPT_LOC .; chmod 755 $SCRIPT_NAME; S3_BUCKET=$S3_BUCKET EXP_IDS=$DOC_IDS NUM_PARTITIONS=$NUM_PA TAXON=$TAXON PARTITION_ID=$p ./$SCRIPT_NAME"
  echo aegea batch submit --command=\"$COMMAND\" --vcpus 16 --memory 64000 --storage /mnt=500 --ecr-image sra_download --environment S3_BUCKET=$S3_BUCKET EXP_IDS=$DOC_IDS NUM_PARTITIONS=$NUM_PA TAXON=$TAXON PARTITION_ID=$p
  echo sleep 20
  p=$(( $p + 1))
done
