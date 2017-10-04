#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Usage: $0 <S3_INPUT_DIR> <S3_OUTPUT_DIR> <TAXON> <CELL_COUNT>"
  echo "   <S3_INPUT_DIR>: S3 Path for the fastq files"
  echo "   <S3_OUTPUT_DIR>: where the results should go"
  echo "   <TAXON>: mus or homo. mus will use mm10 as reference genome, homo will use hg38 as reference genome"
  echo "   <CELL_COUNT>: approx cell count. eg. 3000"
  echo "   <S3_SCRIPT_LOC>: script location i.e. s3://jamestwebber-logs/run_10x_count.py"
  exit 1
fi

S3_INPUT_DIR=$1
S3_OUTPUT_DIR=$2
TAXON=$3
CELL_COUNT=$4
S3_SCRIPT_LOC=$5
SCRIPT_NAME=`basename $S3_SCRIPT_LOC`


COMMAND="aws s3 cp $S3_SCRIPT_LOC .; chmod 755 $SCRIPT_NAME; ./$SCRIPT_NAME"
echo aegea batch submit --command=\"$COMMAND\" --vcpus 32 --memory 128000 --storage /mnt=2000 --ecr-image sra_download --environment S3_INPUT_DIR=$S3_INPUT_DIR S3_OUTPUT_DIR=$S3_OUTPUT_DIR TAXON=$TAXON CELL_COUNT=$CELL_COUNT
