#!/usr/bin/env python
import re
import argparse
import sys
import subprocess
import csv

def generate_sample_map(s3_path, doc_id):
    s3_source = s3_path + '/' + doc_id + '/rawdata/'
    # Generate mapping
    mapping = {}

    # Get sample list
    command = "aws s3 ls %s --recursive | grep fastq.gz" % s3_source
    print command
    output = subprocess.check_output(command, shell=True).split("\n")
    for fname in output:
        m = re.search("/rawdata/([^/]*)/([^/]*)", fname)
        if m:
            sample = m.group(1)
            fastq_name = m.group(2)
            proper_name = fastq_name.split("__")
            mapping[sample] = proper_name[0]

    output_file = "%s.mapping.csv" % doc_id
    with open(output_file, "w") as csvfile:
        m_writer = csv.writer(csvfile, delimiter=',')
        for k,v in mapping.iteritems():
            m_writer.writerow([k, v])
    # copy data to the metadata folder
    command = "aws s3 cp %s %s/%s/metadata/" % (output_file, s3_path, doc_id)
    print command
    output = subprocess.check_output(command, shell=True).split("\n")

    command = "rm -rf %s" % output_file
    print command
    output = subprocess.check_output(command, shell=True).split("\n")


def main():

    description = "Generate mapping from fastq files for Spyros\n"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-s', action="store", dest='s3_path', default=False,
        help="No trailing '/' ")
    parser.add_argument('-d', action="store", dest='exp_id', default=False,
        help="experiment id ")
    results = parser.parse_args()
    if results.s3_path and results.exp_id:
        generate_sample_map(results.s3_path, results.exp_id)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
