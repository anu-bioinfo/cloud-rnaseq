#!/bin/bash
if [ $# -lt 3 ]; then
  echo "Usage: $0 <S3_DIR> <EXP_ID> <SAMPLE_SHEET_FILE_NAME> "
  echo "   <S3_DIR>: S3 Path. Example: s3://czi-hca/data"
  echo "   <EXP_ID>: experiment directory under <S3_DIR> Example: 123456789"
  echo "   <SAMPLE_SHEET_FILE_NAME>: the sample sheet file under <S3_DIR>/<EXP_ID>/config/ Example: sample-sheet.csv"
  exit 1
fi

S3_DIR=$1
EXP_ID=$2
SAMPLE_SHEET=$3

echo aegea batch submit --execute bcl2fastq.py --storage /mnt=2000 --ecr-image sra_download --memory 64000 --environment S3_DIR=${S3_DIR} EXP_ID=${EXP_ID} SAMPLE_SHEET_NAME=${SAMPLE_SHEET} 
