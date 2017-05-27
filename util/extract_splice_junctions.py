#!/usr/bin/env python

# Extracting results/Pass1/SJ.out.tab out

import re
import argparse
import sys
import subprocess
import csv
import os
import threading
import multiprocessing

ROOT_DIR = '/mnt/data'
BATCH_SIZE = 10
def get_doc_list(s3_path):
    command = "aws s3 ls %s/" % s3_path
    print command
    output = subprocess.check_output(command, shell=True)
    doc_list = []
    if output:
        dir_list = output.rstrip().split("\n")
        for f in dir_list:
            m = re.match(".*?([^\/ ]*)/", f)
            if m and len(m.group(1)) > 0:
                doc_list.append(m.group(1))
    return doc_list

def generate_sample_map(s3_source, s3_dest, doc_id, partition_id, num_partitions):
    s3_source = s3_source + '/' + doc_id + '/results/'
    # Generate mapping
    mapping = {}
    generated = []
    count = 0
    # Get sample list
    try:
    #if True:
        command = "aws s3 ls %s --recursive | grep .tgz" % s3_source
        print command
        source_list = subprocess.check_output(command, shell=True).split("\n")

    except Exception:
        print "%s doesn't have the results dir" % doc_id
        return {}

    command = "aws s3 ls %s/%s --recursive |grep .sjout.gz" % (s3_dest, doc_id)
    print command
    try:
        dest_list = subprocess.check_output(command, shell=True).split("\n")
    except Exception:
        dest_list = []
    for f in dest_list:
        m = re.search(doc_id + "\.([^\.]*)\.([^\.]*)\.sjout\.gz", f)
        if m:
            generated.append(m.group(0))
    generated = set(generated)


    for fname in source_list:
        m = re.search("/results/([^\.]*)\.([^\.]*)\.tgz", fname)
        if m:
            if (count % num_partitions) == partition_id:
                sample = m.group(1)
                taxon = m.group(2)
                output_file_name = "%s.%s.%s.sjout.gz" % (doc_id, sample, taxon)
                if output_file_name not in generated:
                    source_file = "%s%s.%s.tgz" % (s3_source, sample, taxon)
                    mapping[output_file_name] = source_file
            count += 1
    return mapping

def run_sample(s3_dest_path, s3_dest_file, s3_source_file):
    m = re.search("/results/([^\.]*\.[^\.]*)\.tgz", s3_source_file)
    if m:
        dir_name = m.group(1)
        data_dir = ROOT_DIR + '/' + dir_name

        command = "mkdir -p %s; aws s3 cp %s %s/" % (data_dir, s3_source_file, data_dir)
        print command
        output = subprocess.check_output(command, shell=True)

        command = "cd %s; tar xvfz %s" % (data_dir, os.path.basename(s3_source_file))
        print command
        output = subprocess.check_output(command, shell=True)

        command = "cd %s; mv results/Pass1/SJ.out.tab %s; gzip %s" % (data_dir, s3_dest_file[:-3], s3_dest_file[:-3])
        print command
        output = subprocess.check_output(command, shell=True)

        command = "cd %s; aws s3 cp %s %s/" % (data_dir, s3_dest_file, s3_dest_path)
        print command
        output = subprocess.check_output(command, shell=True)

        command = "rm -rf %s" % (data_dir)
        print command
        output = subprocess.check_output(command, shell=True)

class sjThread(threading.Thread):
    def __init__(self, run_sample_params):
        threading.Thread.__init__(self)
        self.params = run_sample_params

    def run(self):
        print "%s %s %s" % (self.params[0], self.params[1], self.params[2])
        run_sample(self.params[0], self.params[1], self.params[2])
        sys.stdout.flush()


def run_batch_samples(samples):
    threads = []
    for sample_params in samples:
        sj_thread = sjThread(sample_params)
        sj_thread.start()
        threads.append(sj_thread)
    for t in threads:
        t.join()

def main():

    description = "Harvest splice junction files from the results\n"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-s', action="store", dest='s3_source', default=os.environ.get('S3_SOURCE'))
    parser.add_argument('-o', action="store", dest='s3_dest', default=os.environ.get('S3_DEST'))
    parser.add_argument('-p', action="store", dest='partition', default=os.environ.get('PARTITION'))
    parser.add_argument('-n', action="store", dest='num_partitions', default=os.environ.get('NUM_PARTITIONS'))
    results = parser.parse_args()
    mapping = {}
    if results.s3_source and results.s3_dest and results.partition and results.num_partitions:
        doc_ids = get_doc_list(results.s3_source.rstrip())
        for doc_id in doc_ids:
            doc_map = generate_sample_map(results.s3_source.rstrip(), results.s3_dest.rstrip(),
                    doc_id, int(results.partition), int(results.num_partitions))
            mapping.update(doc_map)
        samples = []
        for k, v in mapping.iteritems():
            print "%s %s" % (k, v)
            samples.append([results.s3_dest, k, v])
            if len(samples) >= BATCH_SIZE:
                run_batch_samples(samples)
                samples = []
        run_batch_samples(samples)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
