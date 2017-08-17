#!/usr/bin/env python
# Example:
#  bcl2fastq --sample-sheet <your sheet> -R <your data dir> -o <fastq_output>
import argparse
import subprocess
import os
import sys
import time

import json
import re

ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + "/data/hca/"
CELLRANGER = 'cellranger'


def main():

    global ROOT_DIR
    global DEST_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        ROOT_DIR = ROOT_DIR + '/' + os.environ['AWS_BATCH_JOB_ID']
        DEST_DIR = ROOT_DIR + '/data/hca/'

    S3_INPUT_DIR  = os.environ['S3_INPUT_DIR'].rstrip('/')
    S3_OUTPUT_DIR  = os.environ['S3_OUTPUT_DIR'].rstrip('/')
    S3_SAMPLE_SHEET_PATH = os.environ['S3_SAMPLE_SHEET_PATH']
    SAMPLE_SHEET_NAME = os.path.basename(S3_SAMPLE_SHEET_PATH)
    EXP_ID = os.path.basename(S3_INPUT_DIR)

    result_path = DEST_DIR + '/' + EXP_ID
    bcl_path = result_path  + '/bcl'
    output_path = result_path + '/fastqs'


    # download sample sheet
    command = "mkdir -p %s; aws s3 cp %s %s/" % (result_path, S3_SAMPLE_SHEET_PATH, result_path)
    print command
    subprocess.check_output(command, shell=True)

    # download the bcl files
    command = "mkdir -p %s; aws s3 cp %s %s --recursive" % (bcl_path, S3_INPUT_DIR, bcl_path)
    print command
    subprocess.check_output(command, shell=True)

    # Run cellranger mkfastq
    command = "%s mkfastq --localmem=48 --sample-sheet=%s/%s --run=%s --output-dir=%s " % (CELLRANGER, result_path, SAMPLE_SHEET_NAME, bcl_path, output_path)
    print command
    subprocess.check_output(command, shell=True)

    # Move fastqs
    command = "aws s3 cp %s %s --recursive" % (output_path, S3_OUTPUT_DIR)
    print command
    subprocess.check_output(command, shell=True)

if __name__ == "__main__":
    main()
