#!/bin/bash

if [ $# -lt 4 ]; then
  echo "Usage: $0 <S3_INPUT_DIR> <S3_OUTPUT_DIR> <TAXON> <CELL_COUNT>"
  echo "   <S3_INPUT_DIR>: S3 Path for the fastq files"
  echo "   <S3_OUTPUT_DIR>: where the results should go"
  echo "   <TAXON>: mus or homo. mus will use mm10 as reference genome, homo will use hg38 as reference genome"
  echo "   <CELL_COUNT>: approx cell count. eg. 3000"
  exit 1
fi

S3_INPUT_DIR=$1
S3_OUTPUT_DIR=$2
TAXON=$3
CELL_COUNT=$4

echo aegea batch submit --execute run_10x_count.py --storage /mnt=500 --ecr-image sra_download --environment S3_INPUT_DIR=$S3_INPUT_DIR S3_OUTPUT_DIR=$S3_OUTPUT_DIR TAXON=$TAXON CELL_COUNT=$CELL_COUNT --memory 64000
