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

def downloaded_to_s3(doc_id):
    s3_dir = S3_BUCKET + '/' + doc_id + '/metadata/'
    try:
        command = "aws s3 ls %s" % (s3_dir)
        #print command
        output = subprocess.check_output(command, shell=True)
        ret_val = []
        if output:
            return True
    except Exception:
        return False


def main():
    global S3_BUCKET
    parser = argparse.ArgumentParser(
        description='Getting list of potential datasets into the HCA dataset')
    parser.add_argument('-s', action="store", dest='s3_path', default=S3_BUCKET)
    parser.add_argument('-f', action="store", dest='sra_doc_list', default=False, help="output file")
    parser.add_argument('-e', action="store", dest='exclude_list', default=False, help='mus or homo')
    parser.add_argument('-o', action="store", dest='output_csv', default = False)
    results = parser.parse_args()
    if results.s3_path and results.sra_doc_list and results.exclude_list and results.output_csv:
        S3_BUCKET=results.s3_path
        with open(results.exclude_list, 'rb') as exclude_file:
            exclude_doc_list = exclude_file.read().split("\n")
        output_list = []
        with open(results.sra_doc_list, 'rb') as csvfile:
           sra_doc_list_reader =  csv.reader(csvfile, delimiter=',')
           for row in sra_doc_list_reader:
               if row[0] in exclude_doc_list:
                   print("%s in exclude_list. ... continue" % row[0])
                   continue
               elif downloaded_to_s3(row[0]):
                   print("%s in is already downloaded. .... continue" % row[0])
                   continue
               # yes a candidate
               print ("%s is a new doc. yay!" % row[0])
               output_list.append(row)

        if len(output_list) > 0:
            # start the mapping
            with open(results.output_csv, 'wb') as csvfile:
                gc_writer = csv.writer(csvfile, delimiter=',')
                for m in output_list:
                    gc_writer.writerow(list(m))
        else:
            print "No new doc_ids for downloading"

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

