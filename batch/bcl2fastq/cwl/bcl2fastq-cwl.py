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
    description = "bcl2fastq wrapper\n"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-d', action="store", dest='s3_path', default=False,
        help="No trailing '/' ")
    parser.add_argument('-s', action="store", dest='sample_sheet_s3_path', default=False,
        help="i.e. s3://biohub-spyros/data/170502/config/sample-sheet.csv")
    results = parser.parse_args()
    if results.s3_path and results.sample_sheet_s3_path:
        s3_source = results.s3_path.rstrip()
        EXP_ID = os.path.basename(s3_source)
        SAMPLE_SHEET_NAME = os.path.basename(results.sample_sheet_s3_path)
    else:
        parser.print_help()
        sys.exit(1)

    global ROOT_DIR
    global DEST_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        ROOT_DIR = ROOT_DIR + '/' + os.environ['AWS_BATCH_JOB_ID']
        DEST_DIR = ROOT_DIR + '/data/hca/'

    result_path = DEST_DIR + EXP_ID
    bcl_path = result_path + '/bcl'
    output_path = result_path + '/fastqs'

    # download sample sheet
    command = "mkdir -p %s; aws s3 cp %s  %s/" % (result_path, results.sample_sheet_s3_path, result_path)
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
