#!/bin/bash
if [ $# -lt 3 ]; then
  echo "Usage: $0 <S3_INPUT_DIR> <S3_OUTPUT_DIR> <S3_SAMPLE_SHEET_PATH>"
  echo "   <S3_INPUT_DIR>: S3 Input Path. Example: s3://czbiohub-seqbot/bcl/170613_A00111_0042_BH27CWDMXX"
  echo "   <S3_OUTPUT_DIR>: S3 Output Path. Example: s3://czbiohub-seqbot/fastqs/170613_A00111_0042_BH27CWDMXX"
  echo "    <S3_SAMPLE_SHEET_PATH>: S3 path for the sample sheet of the bcl. Example: 3://czbiohub-seqbot/sample-sheets/170613_A00111_0042_BH27CWDMXX.csv"
  exit 1
fi

S3_INPUT_DIR=$1
S3_OUTPUT_DIR=$2
S3_SAMPLE_SHEET_PATH=$3

echo aegea batch submit --execute run_10x_mkfastq.py --storage /mnt=2000 --ecr-image sra_download --memory 64000 --environment S3_INPUT_DIR=$S3_INPUT_DIR S3_OUTPUT_DIR=$S3_OUTPUT_DIR S3_SAMPLE_SHEET_PATH=$S3_SAMPLE_SHEET_PATH 
