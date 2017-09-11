#!/bin/bash

# example: ./run_idseq.sh s3://czbiohub-infectious-disease/UGANDA s3://yunfang-workdir/id-uganda UGANDA_S33_NP 999
if [ $# -lt 4 ]; then
  echo "Usage: $0 <INPUT_BUCKET> <OUTPUT_BUCKET> <SAMPLE> <DB_SAMPLE_ID>"
  exit 1
fi

INPUT_BUCKET=$1
OUTPUT_BUCKET=$2
SAMPLE=$3
DB_SAMPLE_ID=$4
S3_SCRIPT_LOC=s3://yunfang-workdir/bin/pipeline.py
SCRIPT_NAME=`basename $S3_SCRIPT_LOC`


COMMAND="aws s3 cp $S3_SCRIPT_LOC .; chmod 755 $SCRIPT_NAME; INPUT_BUCKET=$INPUT_BUCKET OUTPUT_BUCKET=$OUTPUT_BUCKET SAMPLE=$SAMPLE DB_SAMPLE_ID=$DB_SAMPLE_ID ./$SCRIPT_NAME"
  echo aegea batch submit --command=\"$COMMAND\" --storage /mnt=1500 --ecr-image idseq --memory 64000