#!/usr/bin/env python
import json
import time
import os
import sys

from Bio import Entrez
import redis
import redis
import hashlib
import re
from ftplib import FTP
import codecs
import argparse

import GEOparse
import csv
import subprocess

S3_BUCKET = 's3://czi-hca/data'

def get_finish_samples_from_s3(doc_id, taxon):
    s3_dir = S3_BUCKET + '/' + doc_id + '/results/'
    try:
        command = "aws s3 ls %s |grep %s.htseq-count.txt" % (s3_dir, taxon)
        print command
        output = subprocess.check_output(command, shell=True)
        ret_val = []
        if output:
            htseq_files = output.rstrip().split("\n")
            for hf in htseq_files:
                m = re.match(".*?([^\/ ]*)\.(.*)\.htseq-count.txt", hf)
                if m:
                    ret_val.append(m.group(1))
            return ret_val
    except Exception:
        print "%s doesn't have any results" % doc_id
        return []

def get_all_samples_from_s3(doc_id):
    s3_dir = S3_BUCKET + '/' + doc_id + '/rawdata/'
    try:
        command = "aws s3 ls %s" % (s3_dir)
        print command
        output = subprocess.check_output(command, shell=True)
        ret_val = []
        if output:
            sample_dirs = output.rstrip().split("\n")
            for f in sample_dirs:
                m = re.match(".*?([^\/ ]*)/", f)
                if m:
                    ret_val.append(m.group(1))
            return ret_val
    except Exception:
        print "%s doesn't have any samples" % doc_id
        return []


def get_missing_samples_from_s3(doc_id, taxon):
    set_1 = set(get_all_samples_from_s3(doc_id))
    set_2 = set(get_finish_samples_from_s3(doc_id, taxon))
    difference = list(set_1 - set_2)
    ret_val = []
    for sample in difference:
        ret_val.append([doc_id, sample])
    return ret_val

def main():
    global S3_BUCKET
    parser = argparse.ArgumentParser(
        description='Generate gene cell table')
    parser.add_argument('-s', action="store", dest='s3_path', default=False)
    parser.add_argument('-f', action="store", dest='output_csv', default=False)
    parser.add_argument('-t', action="store", dest='taxon', default=False, help='mus or homo')
    parser.add_argument('-d', action="store", dest='doc_ids', default = False)
    results = parser.parse_args()
    if results.s3_path and results.output_csv and results.doc_ids and results.taxon:
        S3_BUCKET=results.s3_path
        missing_list =[]
        doc_list = results.doc_ids.split(',')
        for doc_id in doc_list:
            missing_list += get_missing_samples_from_s3(doc_id, results.taxon)

        if len(missing_list) > 0:
            # start the mapping
            with open(results.output_csv, 'wb') as csvfile:
                gc_writer = csv.writer(csvfile, delimiter=',')
                header_row = ['EXP_ID', 'SAMPLE_ID']
                gc_writer.writerow(header_row)
                for m in missing_list:
                    gc_writer.writerow(list(m))
        else:
            print "No htseq files for %s " % (S3_BUCKET)

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

