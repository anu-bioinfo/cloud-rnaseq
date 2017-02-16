#!/usr/bin/python
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

''' REDIS DB can be downloaded at s3://czsi-sra-config/dataset_info/dump.rdb.gz'''

REDIS_PORT = 7777
REDIS_STORE = redis.StrictRedis(host='localhost', port=REDIS_PORT, db=2)
FTP_HOST = "ftp://ftp-trace.ncbi.nlm.nih.gov"
DEST_DIR = "/mnt/data/hca/"
WGET = '/usr/bin/wget '
MIN_FILE_SIZE = 5000
S3_BUCKET = 's3://czi-hca/data'


STAR="/usr/local/bin/STAR"
HTSEQ="/usr/local/bin/htseq-count"
SAMTOOLS="/usr/local/bin/samtools"

GENOME_DIR="/mnt/genome/STAR/HG38/" # change
GENOME_FASTA="/mnt/genome/hg38/hg38.fa" # change
SJDB_GTF="/mnt/genome/hg38/hg38.gtf" # change

COMMON_PARS="--runThreadN 12 --outFilterType BySJout \
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
--readFilesCommand zcat"


def run_sample(sample_name, doc_id):
	''' Example:
	   run_sample(SRR1974579, 200067835)
	'''
	dest_dir = DEST_DIR + sample_name
	subprocess.check_output("mkdir -p %s" % dest_dir, shell=True)
	subprocess.check_output("mkdir -p %s/rawdata" % dest_dir, shell=True)
	subprocess.check_output("mkdir -p %s/results" % dest_dir, shell=True)

	# copy fast.gz from s3 to local
	s3_source = S3_BUCKET + '/' + doc_id + '/rawdata/' + sample_name + '/'
	command = "aws s3 cp %s %s/rawdata/ --recursive --recursive --exclude '*' --include '*.fastq.gz'" % (s3_source, dest_dir) 
	print command
	output = subprocess.check_output(command, shell=True)

	# start running STAR
	# getting input files first
	command = "ls %s/rawdata/*.fastq.gz" % dest_dir
	print command
	reads = subprocess.check_output(command, shell=True).replace("\n", " ")
	if len(reads) < 20:
		print "Empty reads for %s" % s3_source
		return
	command = "cd %s/results; mkdir -p Pass1; cd Pass1;" % dest_dir
	command += STAR + ' ' + COMMON_PARS + ' --genomeDir ' + GENOME_DIR + ' --sjdbGTFfile ' + SJDB_GTF + ' --readFilesIn ' + reads
	print command
	output = subprocess.check_output(command, shell=True)
	# running sam tools
	command = "cd %s/results; %s sort -n -m 6000000000 -o Pass1/Aligned.out.sorted.bam Pass1/Aligned.out.bam" % (dest_dir, SAMTOOLS)
	print command
	output = subprocess.check_output(command, shell=True)

	# running htseq
	command = "cd %s/results; %s -s no -f bam -m intersection-nonempty  ./Pass1/Aligned.out.sorted.bam  %s > htseq-count.txt" % (dest_dir, HTSEQ, SJDB_GTF)
	print command
	output = subprocess.check_output(command, shell=True)

	# compressed the results dir and move it to s3
	command = "cd %s; tar cvfz %s.tgz results" % (dest_dir, sample_name)
	print command
	output = subprocess.check_output(command, shell=True)

	# copy htseq and log files out to s3
	s3_dest = S3_BUCKET + '/' + doc_id + '/results/'
	command = "aws s3 cp %s/%s.tgz %s" % (dest_dir, sample_name, s3_dest)
	print command
	subprocess.check_output(command, shell=True)

	command = "aws s3 cp %s/results/htseq-count.txt %s%s.htseq-count.txt" % (dest_dir, s3_dest, sample_name)
	print command
	subprocess.check_output(command, shell=True)

	command = "aws s3 cp %s/results/Pass1/Log.final.out %s%s.log.final.out" % (dest_dir, s3_dest, sample_name)
	print command
	subprocess.check_output(command, shell=True)

	# rm all the files
	command = "rm -rf %s" % dest_dir
	print command
	subprocess.check_output(command, shell=True)
	
	

def run(doc_id, num_partitions, partition_id):
	print "Running partition %d of %d for doc %s" % (partition_id, num_partitions, doc_id)
	command = "aws s3 ls %s/%s/rawdata/" % (S3_BUCKET, doc_id)
	print command
	output = subprocess.check_output(command, shell=True).split("\n")
	sample_list = []
	
	for f in output:
		matched = re.search("\s([\d\w]+)\/",f)
		if matched: 
			sample_list.append(matched.group(1))
	idx = 0
	for sample_name in sample_list:
		if idx % num_partitions == partition_id:
			try:
				run_sample(sample_name, doc_id)
				print "%s : %s " % (doc_id, sample_name)
			except Exception:
				print "Error downloading %s. Retry in 5 seconds" % sample_name
				
		idx += 1


def main():
	# copy the redis database and start the redis-server
	doc_id = os.environ['GDS_ID']
	num_partitions = int(os.environ['NUM_PARTITIONS'])
	partition_id = int(os.environ['PARTITION_ID'])

	# download the genome data
	command = "mkdir -p /mnt/genome; aws s3 cp s3://czi-hca/ref-genome/hg38.tgz /mnt/genome"
	print command
	subprocess.check_output(command, shell=True)

	command = "cd /mnt/genome; tar xvfz hg38.tgz"
	print command
	subprocess.check_output(command, shell=True)

	command = "mkdir -p /mnt/genome/STAR; aws s3 cp s3://czi-hca/ref-genome/STAR/HG38.tgz /mnt/genome/STAR"
	print command
	subprocess.check_output(command, shell=True)

	command = "cd /mnt/genome/STAR; tar xvfz HG38.tgz"
	print command
	subprocess.check_output(command, shell=True)

	# run through the samples
	run(doc_id, num_partitions, partition_id)


if __name__ == "__main__":
	main()
