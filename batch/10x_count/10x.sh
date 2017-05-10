#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: $0 <S3_DIR> <TAXON> <CELL_COUNT>"
  echo "   <S3_DIR>: S3 Path for the library. Example: s3://czi-hca/data"
  echo "   <TAXON>: mus or homo. mus will use mm10 as reference genome, homo will use hg38 as reference genome"
  echo "   <CELL_COUNT>: approx cell count. eg. 3000"
  exit 1
fi

S3_DIR=$1
TAXON=$2
CELL_COUNT=$3

echo aegea batch submit --execute run_10x_count.py --storage /mnt=500 --ecr-image sra_download --environment S3_DIR=$S3_DIR TAXON=$TAXON CELL_COUNT=$CELL_COUNT --memory 64000
