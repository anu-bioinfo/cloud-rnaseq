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

CELL_RANGER = 'cellranger'
GC_TABLE_GENERATOR = 'generate_gc_table_from_cellranger.py'

ROOT_DIR = '/mnt'


def main():
    root_dir = ROOT_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        root_dir = os.path.join(root_dir, os.environ['AWS_BATCH_JOB_ID'])

    dest_dir = os.path.join(root_dir, 'data', 'hca')

    # required environmental variable
    taxon = os.environ.get('TAXON', 'homo')
    s3_input_dir = os.environ['S3_INPUT_DIR'].rstrip('/')
    s3_output_dir = os.environ['S3_OUTPUT_DIR'].rstrip('/')
    cell_count = os.environ['CELL_COUNT']


    genome_base_dir = os.path.join(root_dir, "genome/cellranger")
    if taxon == 'homo':
        genome_name = 'HG38-PLUS'
    elif taxon == 'mus':
        genome_name = 'MM10-PLUS'
    else:
        raise ValueError("unknown taxon {}".format(taxon))

    # files that should be uploaded outside of the massive tgz
    # path should be relative to the run folder
    files_to_upload = [
        'outs/raw_gene_bc_matrices_h5.h5',
        'outs/raw_gene_bc_matrices/{}/genes.tsv'.format(genome_name),
        'outs/raw_gene_bc_matrices/{}/barcodes.tsv'.format(genome_name),
        'outs/raw_gene_bc_matrices/{}/matrix.mtx'.format(genome_name),
        'outs/web_summary.html',
        'outs/metrics_summary.csv'
    ]

    genome_tar_source = os.path.join('s3://czi-hca/ref-genome/cellranger/',
                                     genome_name + '.tgz')
    genome_dir = os.path.join(genome_base_dir, genome_name)

    if os.environ.get('AWS_BATCH_JOB_ID'): # Using BATCH
        # download the ref genome data
        os.makedirs(genome_base_dir)
        command = "aws s3 cp {} {}/".format(genome_tar_source, genome_base_dir)
        print command
        subprocess.check_output(command, shell=True)

        ref_genome_file = os.path.basename(genome_tar_source)
        os.chdir(genome_base_dir)
        command = "tar xvfz {}".format(ref_genome_file)
        print command
        subprocess.check_output(command, shell=True)

    sys.stdout.flush()

    # download the fastq files
    sample_id = os.path.basename(s3_input_dir)
    result_path = os.path.join(dest_dir, sample_id)
    fastq_path = os.path.join(result_path, 'fastqs')


    os.makedirs(fastq_path)
    command = "aws s3 cp {} {} --recursive".format(s3_input_dir, fastq_path)
    print command
    subprocess.check_output(command, shell=True)


    # Run cellranger
    os.chdir(result_path)
    command = ("{} count --localmem=40 --nosecondary --cells={}"
               " --sample={} --id={} --fastqs={} --transcriptome={}").format(
            CELL_RANGER, cell_count, sample_id,
            sample_id, fastq_path, genome_dir
    )
    print command
    try:
        cmnd_output = subprocess.check_output(command,
                                              stderr=subprocess.STDOUT,
                                              shell=True,
                                              universal_newlines=True)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        command = "tar cvfz {}.tgz {}".format(sample_id, sample_id)
        print command
        subprocess.check_output(command, shell=True)

        command = "aws s3 cp {}.tgz {}/".format(os.path.join(result_path,
                                                             sample_id),
                                                s3_output_dir)
        print command
        subprocess.check_output(command, shell=True)
        sys.exit(1)

    # Move results(websummary, cell-gene table, tarball) data back to S3
    for file_name in files_to_upload:
        command = "aws s3 cp {} {}/".format(
                os.path.join(result_path, sample_id, file_name), s3_output_dir
        )
        print command
        subprocess.check_output(command, shell=True)


    command = "tar cvfz {}.tgz {}".format(
            os.path.join(result_path, sample_id), sample_id
    )
    print command
    subprocess.check_output(command, shell=True)


    command = "aws s3 cp {}.tgz {}/ ".format(
            os.path.join(result_path, sample_id), s3_output_dir
    )
    print command
    subprocess.check_output(command, shell=True)


    matrix_dir = "raw_gene_bc_matrices/" + genome_name
    command = "{} -d {} -f {}.{}.cell-gene.csv -m 500".format(
            GC_TABLE_GENERATOR, os.path.join(sample_id, 'outs', matrix_dir),
            sample_id, taxon
    )
    print command
    subprocess.check_output(command, shell=True)


    command = "aws s3 cp {}.{}.cell-gene.csv {}/".format(
            sample_id, taxon, s3_output_dir
    )
    print command
    subprocess.check_output(command, shell=True)


if __name__ == "__main__":
    main()
