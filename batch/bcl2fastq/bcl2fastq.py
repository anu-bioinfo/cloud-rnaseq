#!/usr/bin/env python
# Example:
#  bcl2fastq --sample-sheet <your sheet> -R <your data dir> -o <fastq_output>
import argparse
import subprocess
import os
import sys
import thread
import time

from Bio import Entrez

import json
import redis
import re

import threading
import multiprocessing

ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + "/data/hca/"
BCL2FASTQ = 'bcl2fastq'

def main():

    global ROOT_DIR
    global DEST_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        ROOT_DIR = ROOT_DIR + '/' + os.environ['AWS_BATCH_JOB_ID']
        DEST_DIR = ROOT_DIR + '/data/hca/'

    S3_DIR  = os.environ['S3_DIR'].rstrip('/')
    EXP_ID = os.environ['EXP_ID'].split(",")
    SAMPLE_SHEET_NAME = os.environ.get('SAMPLE_SHEET', 'sample-sheet.csv')

    result_path = DEST_DIR + EXP_ID
    bcl_path = result_path + '/bcl'
    output_path = result_path + '/fastqs'
    s3_source = S3_DIR + '/' + EXP_ID

    # download sample sheet
    command = "mkdir -p %s; aws s3 cp %s/config/%s %s/" % (result_path, s3_source, SAMPLE_SHEET_NAME, result_path)
    print command
    subprocess.check_output(command, shell=True)

    # download the bcl files
    command = "mkdir -p %s; aws s3 cp %s/bcl/ %s --recursive" % (bcl_path, s3_source, bcl_path)
    print command
    subprocess.check_output(command, shell=True)

    # Run bcl2 fastq
    command = "%s --sample-sheet %s/%s -R %s -o %s " % (BCL2FASTQ, result_path, SAMPLE_SHEET_NAME, bcl_path, output_path)
    print command
    subprocess.check_output(command, shell=True)

    # Move reports(websummary, cell-gene table, tarball) data back to S3
    reports_path = subprocess.check_output("ls -d %s/Reports/html/*/all/all/all/" % output_path, shell=True).rstrip()
    command = "aws s3 cp %s %s/reports/bcl2fastq/ --recursive " % (reports_path, s3_source)
    print command
    subprocess.check_output(command, shell=True)

    # Move fastqs
    fastq_dir_command = "ls %s |grep -v '.gz' |grep -v Reports|grep -v Stats" % output_path
    fastq_dir_names = subprocess.check_output(fastq_dir_command, shell = True).rstrip().split("\n")

    for fastq_dir in  fastq_dir_names:
        command = "aws s3 cp %s/%s/ %s/rawdata/ --recursive" % (output_path, fastq_dir, s3_source)
        print command
        subprocess.check_output(command, shell=True)

if __name__ == "__main__":
    main()
