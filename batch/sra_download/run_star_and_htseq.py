#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import thread
import time
import datetime

from Bio import Entrez

import json
import redis
import re

import threading
import multiprocessing

import logging

ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + "/data/hca/"
S3_BUCKET = 's3://czi-hca/data'


STAR="/usr/local/bin/STAR"
HTSEQ="htseq-count"
SAMTOOLS="samtools"

GENOME_DIR = ROOT_DIR + "/genome/STAR/HG38-PLUS/" # change
GENOME_FASTA = ROOT_DIR + "/genome/hg38-plus/hg38-plus.fa" # change
SJDB_GTF = ROOT_DIR + "/genome/hg38-plus/hg38-plus.gtf" # change

TAXON="homo"

COMMON_PARS="--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI NM MD \
--genomeLoad LoadAndKeep \
--outReadsUnmapped Fastx \
--readFilesCommand zcat"

HTSEQ_THREADS_MAX = 4
CURR_MIN_VER = datetime.date(2017, 3, 1)

def maybe_log(msg, logger=None, **args):
    if logger is not None:
        logger.info(msg, **args)
    print msg


def run_sample(sample_name, doc_id, logger=None):
    """Example:
       run_sample(SRR1974579, 200067835)
    """

    dest_dir = DEST_DIR + sample_name
    subprocess.check_output("mkdir -p %s" % dest_dir, shell=True)
    subprocess.check_output("mkdir -p %s/rawdata" % dest_dir, shell=True)
    subprocess.check_output("mkdir -p %s/results" % dest_dir, shell=True)

    # copy fast.gz from s3 to local
    s3_source = S3_BUCKET + '/' + doc_id + '/rawdata/' + sample_name + '/'
    command = "aws s3 cp %s %s/rawdata/ --recursive --exclude '*' --include '*.fastq.gz'" % (s3_source, dest_dir)
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)

    # start running STAR
    # getting input files first
    command = "ls %s/rawdata/*.fastq.gz" % dest_dir
    maybe_log(command, logger)
    reads = subprocess.check_output(command, shell=True).replace("\n", " ")
    if len(reads) < 20:
        maybe_log("Empty reads for %s" % s3_source)
        return

    command = "cd %s/results; mkdir -p Pass1; cd Pass1;" % dest_dir
    #command += STAR + ' ' + COMMON_PARS + ' --genomeDir ' + GENOME_DIR + ' --sjdbGTFfile ' + SJDB_GTF + ' --readFilesIn ' + reads
    command += STAR + ' ' + COMMON_PARS + ' --runThreadN ' + str(multiprocessing.cpu_count()) +  ' --genomeDir ' + GENOME_DIR + ' --readFilesIn ' + reads
    maybe_log(command, logger)
    sys.stdout.flush()
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)

    # running sam tools
    command = "cd %s/results; %s sort -m 6000000000 -o ./Pass1/Aligned.out.sorted.bam ./Pass1/Aligned.out.bam" % (dest_dir, SAMTOOLS)
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)
    sys.stdout.flush()

    # running samtools index -b
    command = "cd %s/results/Pass1; %s index -b Aligned.out.sorted.bam " % (dest_dir, SAMTOOLS)
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)
    sys.stdout.flush()

    # remove unsorted bam files
    command = "cd %s/results/Pass1; rm -rf Aligned.out.bam " % dest_dir
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)


    # remove fastq files
    command = "rm -rf %s/rawdata/*.fastq.gz" % dest_dir
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)

    # generating files for htseq-count
    command = "cd %s/results; %s sort -m 6000000000 -n -o ./Pass1/Aligned.out.sorted-byname.bam  ./Pass1/Aligned.out.sorted.bam " % (dest_dir, SAMTOOLS)
    maybe_log(command, logger)
    output = subprocess.check_output(command, shell=True)
    maybe_log(output, logger)

    sys.stdout.flush()

    # ready to be htseq-ed and cleaned up
    return { 'doc_id': doc_id, 'sample_name': sample_name, 'dest_dir': dest_dir}


class htseqThread(threading.Thread):
    def __init__(self, htseqParams, logger=None):
        threading.Thread.__init__(self)
        self.htseqParams = htseqParams
        self.logger = logger

    def run(self):
        doc_id = self.htseqParams['doc_id']
        sample_name = self.htseqParams['sample_name']
        dest_dir = self.htseqParams['dest_dir']

        # running htseq

        command = "cd %s/results; %s -r name -s no -f bam -m intersection-nonempty  ./Pass1/Aligned.out.sorted-byname.bam  %s > htseq-count.txt; rm ./Pass1/Aligned.out.sorted-byname.bam" % (dest_dir, HTSEQ, SJDB_GTF)
        maybe_log(command, self.logger)
        output = subprocess.check_output(command, shell=True)
        maybe_log(output, self.logger)

        # compressed the results dir and move it to s3
        command = "cd %s; tar cvfz %s.%s.tgz results" % (dest_dir, sample_name, TAXON)
        maybe_log(command, self.logger)
        output = subprocess.check_output(command, shell=True)
        maybe_log(output, self.logger)

        # copy htseq and log files out to s3
        s3_dest = S3_BUCKET + '/' + doc_id + '/results/'
        command = "aws s3 cp %s/%s.%s.tgz %s" % (dest_dir, sample_name, TAXON, s3_dest)
        maybe_log(command, self.logger)
        output = subprocess.check_output(command, shell=True)
        maybe_log(output, self.logger)

        command = "aws s3 cp %s/results/htseq-count.txt %s%s.%s.htseq-count.txt" % (dest_dir, s3_dest, sample_name, TAXON)
        maybe_log(command, self.logger)
        output = subprocess.check_output(command, shell=True)
        maybe_log(output, self.logger)

        command = "aws s3 cp %s/results/Pass1/Log.final.out %s%s.%s.log.final.out" % (dest_dir, s3_dest, sample_name, TAXON)
        maybe_log(command, self.logger)
        output = subprocess.check_output(command, shell=True)
        maybe_log(output, self.logger)

        # rm all the files
        command = "rm -rf %s" % dest_dir
        maybe_log(command, self.logger)
        subprocess.check_output(command, shell=True)

        sys.stdout.flush()


def runHtseq(htseq_jobs, logger):
    threads = []

    for htseqParams in htseq_jobs:
        maybe_log("Starting a htseq thread for %s " % htseqParams, logger)
        ht_thread = htseqThread(htseqParams, logger)
        ht_thread.start()
        threads.append(ht_thread)

    for t in threads:
        t.join()


def run(doc_ids, num_partitions, partition_id, logger=None):
    htseq_jobs = []
    for doc_id in doc_ids:
        if doc_id.startswith('Undetermined'):
            maybe_log("Skipping file: %s" % doc_id, logger)
            continue

        # Check the doc_id folder for existing runs
        command = "aws s3 ls {}/".format(
                os.path.join(S3_BUCKET, doc_id, 'results')
        )
        maybe_log(command, logger)
        output = subprocess.check_output(command, shell=True)

        output_files = {(line[:10].split('-'), line.split()[-1])
                        for line in output.split('\n')
                        if line.strip().endswith('htseq-count.txt')}
        output_files = {fn for dt,fn in output_files
                        if datetime.date(*map(int, dt)) > CURR_MIN_VER}
        maybe_log("number of files: {}".format(len(output_files)), logger)

        maybe_log("Running partition {} of {} for doc {}".format(
                partition_id, num_partitions, doc_id
        ), logger)
        command = "aws s3 ls %s/%s/rawdata/" % (S3_BUCKET, doc_id)
        maybe_log(command, logger)
        try:
            output = subprocess.check_output(command, shell=True).split("\n")
        except subprocess.CalledProcessError:
            maybe_log("error in aws command", logger, exc_info=True)
            output = []

        sample_list = []

        for f in output:
            matched = re.search("\s([\d\w\-.]+)/", f)
            if matched:
                sample_list.append(matched.group(1))

        for sample_name in sample_list[partition_id::num_partitions]:
            if '{}.{}.htseq-count.txt'.format(sample_name, TAXON) in output_files:
                maybe_log("{} already exists, skipping".format(sample_name, logger))
                continue

            try:
                maybe_log("Running sample {}".format(sample_name), logger)
                ret = run_sample(sample_name, doc_id)
                maybe_log("%s : %s " % (doc_id, sample_name), logger)
                if ret is not None:
                    htseq_jobs.append(ret)
                if len(htseq_jobs) >= HTSEQ_THREADS_MAX:
                    runHtseq(htseq_jobs, logger)
                    htseq_jobs = []

            except subprocess.CalledProcessError, e:
                maybe_log("Error processing stdout output: %s for sample %s\n" % (e.output, sample_name), logger)
            except:
                maybe_log("Error processing %s. " % sample_name, logger, exc_info=True)

    if htseq_jobs:
        runHtseq(htseq_jobs, logger)


def main(logger=None):
    global ROOT_DIR
    global DEST_DIR

    if os.environ.get('AWS_BATCH_JOB_ID'):
        ROOT_DIR = ROOT_DIR + '/' + os.environ['AWS_BATCH_JOB_ID']
        DEST_DIR = ROOT_DIR + '/data/hca/'

    global GENOME_DIR
    global GENOME_FASTA
    global SJDB_GTF
    global TAXON

    taxon = os.environ.get('TAXON', 'homo')
    if taxon == 'homo':
        GENOME_DIR = ROOT_DIR + "/genome/STAR/HG38-PLUS/" # change
        GENOME_FASTA = ROOT_DIR + "/genome/hg38-plus/hg38-plus.fa" # change
        SJDB_GTF = ROOT_DIR + "/genome/hg38-plus/hg38-plus.gtf" # change
        TAXON = 'homo'

        ref_genome_file = 'hg38-plus.tgz'
        ref_genome_star_file = 'HG38-PLUS.tgz'
    elif taxon == 'mus':
        GENOME_DIR = ROOT_DIR + "/genome/STAR/MM10-PLUS/"
        GENOME_FASTA = ROOT_DIR + "/genome/mm10-plus/mm10-plus.fa" # change (?)
        SJDB_GTF = ROOT_DIR + "/genome/mm10-plus/mm10-plus.gtf" # change (?)
        TAXON = 'mus'

        ref_genome_file = 'mm10-plus.tgz'
        ref_genome_star_file = 'MM10-PLUS.tgz'
    else:
        raise ValueError('Invalid taxon {}'.format(taxon))

    maybe_log(
            'Run Info:\n{}\n'.format(
                    '\n'.join((GENOME_DIR, GENOME_FASTA, SJDB_GTF, TAXON,
                               ref_genome_file, ref_genome_star_file,
                               os.environ['S3_BUCKET'],
                               os.environ['EXP_IDS'],
                               os.environ['NUM_PARTITIONS'],
                               os.environ['PARTITION_ID']))
            ),
            logger
    )

    # download the genome data

    command = "mkdir -p %s/genome; aws s3 cp s3://czi-hca/ref-genome/%s %s/genome" % (ROOT_DIR, ref_genome_file, ROOT_DIR)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)

    command = "cd %s/genome; tar xvfz %s " % (ROOT_DIR, ref_genome_file)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)

    command = "mkdir -p %s/genome/STAR; aws s3 cp s3://czi-hca/ref-genome/STAR/%s %s/genome/STAR" % (ROOT_DIR, ref_genome_star_file, ROOT_DIR)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)

    command = "cd %s/genome/STAR; tar xvfz %s" % (ROOT_DIR, ref_genome_star_file)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)


    sys.stdout.flush()

    # Load Genome Into Memory
    command = "%s --genomeDir %s --genomeLoad LoadAndExit" % (STAR, GENOME_DIR)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)

    # run through the samples
    global S3_BUCKET
    S3_BUCKET = os.environ['S3_BUCKET']

    num_partitions = int(os.environ['NUM_PARTITIONS'])
    partition_id = int(os.environ['PARTITION_ID'])
    doc_ids = os.environ['EXP_IDS'].split(",")
    run(doc_ids, num_partitions, partition_id, logger)

    # Remove Genome from Memory
    command = "%s --genomeDir %s --genomeLoad Remove" % (STAR, GENOME_DIR)
    maybe_log(command, logger)
    subprocess.check_output(command, shell=True)


if __name__ == "__main__":
    if os.environ.get('AWS_BATCH_JOB_ID'):
        logging.basicConfig(level=logging.INFO)
        mainlogger = logging.getLogger(__name__)

        log_file = os.path.join(
                ROOT_DIR, '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
        )
        handler = logging.FileHandler(log_file)
        handler.setLevel(mainlogger.level)

        # create a logging format
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)

        # add the handlers to the logger
        mainlogger.addHandler(handler)

        try:
            main(mainlogger)
        finally:
            command = 'aws s3 cp {} s3://jamestwebber-logs/star_logs/'.format(log_file)
            mainlogger.info(command)
            print command

            handler.close()
            subprocess.check_output(command, shell=True)
    else:
        main()
