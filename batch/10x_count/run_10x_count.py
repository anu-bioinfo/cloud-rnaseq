#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./run_10x_count.py
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

#CELL_RANGER = '/cellranger-1.3.1/cellranger'
CELL_RANGER = 'cellranger'
GC_TABLE_GENERATOR = 'generate_gc_table_from_cellranger.py'

ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + "/data/hca/"

def main():

    global ROOT_DIR
    global DEST_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        ROOT_DIR = ROOT_DIR + '/' + os.environ['AWS_BATCH_JOB_ID']
        DEST_DIR = ROOT_DIR + '/data/hca/'

    # required environmental variable
    TAXON   = os.environ['TAXON']
    S3_INPUT_DIR  = os.environ['S3_INPUT_DIR'].rstrip('/')
    S3_OUTPUT_DIR = os.environ['S3_OUTPUT_DIR'].rstrip('/')
    CELL_COUNT = os.environ['CELL_COUNT']


    GENOME_BASE_DIR = ROOT_DIR + "/genome/cellranger"
    GENOME_NAME = 'HG38-PLUS'
    if TAXON == 'mus':
        GENOME_NAME = 'MM10-PLUS'

    GENOME_TAR_SORUCE = 's3://czi-hca/ref-genome/cellranger/' + GENOME_NAME + '.tgz'
    GENOME_DIR = GENOME_BASE_DIR + "/" + GENOME_NAME + "/"

    # download the ref genome data
    command = "mkdir -p %s; aws s3 cp %s  %s/" %  (GENOME_BASE_DIR, GENOME_TAR_SORUCE, GENOME_BASE_DIR)
    print command
    #subprocess.check_output(command, shell=True)

    ref_genome_file = os.path.basename(GENOME_TAR_SORUCE)
    command = "cd %s; tar xvfz %s" % (GENOME_BASE_DIR, ref_genome_file)
    print command
    #subprocess.check_output(command, shell=True)

    sys.stdout.flush()

    # download the fastq files
    sample_id = os.path.basename(S3_INPUT_DIR)
    result_path = DEST_DIR + sample_id
    fastq_path = result_path + '/fastqs'
    output_path = result_path + '/' + sample_id


    command = "mkdir -p %s; aws s3 cp %s %s --recursive" % (fastq_path, S3_INPUT_DIR, fastq_path)
    print command
    subprocess.check_output(command, shell=True)

    # Run cellranger
    command = "cd %s; %s count --cells=%s --sample=%s --id=%s --fastqs=%s --transcriptome=%s" % (result_path, CELL_RANGER, CELL_COUNT, sample_id, sample_id, fastq_path, GENOME_DIR)
    print command
    subprocess.check_output(command, shell=True)

    # Move results(websummary, cell-gene table, tarball) data back to S3
    command = "cd %s; aws s3 cp %s/outs/web_summary.html %s/ " % (result_path, sample_id, S3_OUTPUT_DIR)
    print command
    subprocess.check_output(command, shell=True)


    command = "cd %s; tar cvfz %s.tgz %s" % (result_path, sample_id, sample_id)
    print command
    subprocess.check_output(command, shell=True)

    command = "cd %s; aws s3 cp %s.tgz %s/ " % (result_path, sample_id, S3_OUTPUT_DIR)
    print command
    subprocess.check_output(command, shell=True)

    matrix_dir = "raw_gene_bc_matrices/" + GENOME_NAME
    command = "cd %s; %s -d %s/outs/%s -f %s.%s.cell-gene.csv -m 500" % (result_path, GC_TABLE_GENERATOR, sample_id, matrix_dir, sample_id, TAXON)
    print command
    subprocess.check_output(command, shell=True)

    command = "cd %s; aws s3 cp %s.%s.cell-gene.csv %s/" % (result_path, sample_id, TAXON, S3_OUTPUT_DIR)
    print command
    subprocess.check_output(command, shell=True)

if __name__ == "__main__":
    main()
