#!/usr/bin/env python
import re
import argparse
import sys
import subprocess
import csv

def move_fastqs_to_dirs(s3_path, doc_id, s3_outputpath):
    s3_source = s3_path.rstrip('/') + '/' + doc_id # + '/rawdata'
    s3_output_source = s3_outputpath.rstrip('/') + '/' + doc_id + '/rawdata'
    # Generate mapping
    mapping = {}

    # Get sample list
    command = "aws s3 ls %s --recursive | grep fastq.gz" % s3_source
    print command
    output = subprocess.check_output(command, shell=True).split("\n")
    for fname in output:
        m = re.search("/([^/]*)(_R[12]_001.fastq.gz)", fname)
        if m:
            sample = m.group(1)
            fastq_name = m.group(1) + m.group(2)
            command = "aws s3 cp %s/%s %s/%s/ " % (s3_source, fastq_name, s3_output_source, sample)
            print command
            output = subprocess.check_output(command, shell=True).split("\n")


def main():

    description = "Move Fastqs from flat to directory structures\n"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-s', action="store", dest='s3_path', default=False,
        help="No trailing '/' ")
    parser.add_argument('-d', action="store", dest='exp_id', default=False,
        help="experiment id ")
    parser.add_argument('-o', action="store", dest='s3_output_path', default=False,
        help="output s3 path")
    results = parser.parse_args()
    if results.s3_path and results.exp_id and results.s3_output_path:
        move_fastqs_to_dirs(results.s3_path, results.exp_id, results.s3_output_path)
    else:
        print "Example:./fastq_move.py -s s3://czbiohub-seqbot/fastqs -d 170824_NB501961_0009_AHMJFTBGX2 -o s3://czbiohub-maca/data"
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
