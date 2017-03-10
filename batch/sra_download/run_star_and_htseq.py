#!/usr/bin/env python
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


DEST_DIR = "/mnt/data/hca/"
S3_BUCKET = 's3://czi-hca/data'


STAR="/usr/local/bin/STAR"
HTSEQ="htseq-count"
SAMTOOLS="samtools"

GENOME_DIR="/mnt/genome/STAR/HG38-PLUS/" # change
GENOME_FASTA="/mnt/genome/hg38-plus/hg38-plus.fa" # change
SJDB_GTF="/mnt/genome/hg38-plus/hg38-plus.gtf" # change

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
--genomeLoad LoadAndKeep \
--outReadsUnmapped Fastx \
--readFilesCommand zcat"

HTSEQ_THREADS_MAX = 6
CURR_MIN_VER = 20170301

def run_sample(sample_name, doc_id, force_download = False):
	''' Example:
	   run_sample(SRR1974579, 200067835)
	'''


	# Check if star has been run
	if not force_download:
		s3_dest = S3_BUCKET + '/' + doc_id + '/results/'
		command = "aws s3 ls %s | grep %s.htseq" % (s3_dest, sample_name)
                print command
		try:
			output = subprocess.check_output(command, shell=True)
			if len(output) > 10:
				# sample has been downloaded before
				dateint = int(output[0:10].replace('-', ''))
				if dateint > CURR_MIN_VER:
					return
		except Exception:
			# No big deal. S3 error
			print "%s doesn't exist yet" % sample_name

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
	#command += STAR + ' ' + COMMON_PARS + ' --genomeDir ' + GENOME_DIR + ' --sjdbGTFfile ' + SJDB_GTF + ' --readFilesIn ' + reads
	command += STAR + ' ' + COMMON_PARS + ' --genomeDir ' + GENOME_DIR + ' --readFilesIn ' + reads
	print command
	sys.stdout.flush()
	output = subprocess.check_output(command, shell=True)
	# running sam tools
	command = "cd %s/results; %s sort -m 6000000000 -o ./Pass1/Aligned.out.sorted.bam ./Pass1/Aligned.out.bam" % (dest_dir, SAMTOOLS)
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	sys.stdout.flush()

	# running samtools index -b
	command = "cd %s/results/Pass1; %s index -b Aligned.out.sorted.bam " % (dest_dir, SAMTOOLS)
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	sys.stdout.flush()

	# remove unsorted bam files
	command = "cd %s/results/Pass1; rm -rf Aligned.out.bam " % (dest_dir)
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	sys.stdout.flush()

	# ready to be htseq-ed and cleaned up
	return { 'doc_id': doc_id, 'sample_name': sample_name, 'dest_dir': dest_dir}


class htseqThread(threading.Thread):
	def __init__(self, htseqParams):
		threading.Thread.__init__(self)
		self.htseqParams = htseqParams

	def run(self):
		doc_id = self.htseqParams['doc_id']
		sample_name = self.htseqParams['sample_name']
		dest_dir = self.htseqParams['dest_dir']

		# running htseq
		command = "cd %s/results; %s -r pos -s no -f bam -m intersection-nonempty  ./Pass1/Aligned.out.sorted.bam  %s > htseq-count.txt" % (dest_dir, HTSEQ, SJDB_GTF)
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

		sys.stdout.flush()

def runHtseq(htseq_jobs):
	threads = []
	for htseqParams in htseq_jobs:
		print "Starting a htseq thread for %s " % htseqParams
		ht_thread = htseqThread(htseqParams)
		ht_thread.start()
		threads.append(ht_thread)
	for t in threads:
		t.join()

def run(doc_ids	, num_partitions, partition_id):
	htseq_jobs = []
	for doc_id in doc_ids:
		print "Running partition %d of %d for doc %s" % (partition_id, num_partitions, doc_id)
		command = "aws s3 ls %s/%s/rawdata/" % (S3_BUCKET, doc_id)
		print command
		output = subprocess.check_output(command, shell=True).split("\n")
		sample_list = []

		for f in output:
			matched = re.search("\s([\d\w\-\.]+)\/",f)
			if matched:
				sample_list.append(matched.group(1))
		idx = 0
		for sample_name in sample_list:
			if idx % num_partitions == partition_id:
				try:
				#if True:
					ret = run_sample(sample_name, doc_id)
					print "%s : %s " % (doc_id, sample_name)
					if ret is not None:
						htseq_jobs.append(ret)
					if len(htseq_jobs) >= HTSEQ_THREADS_MAX:
						runHtseq(htseq_jobs)
						htseq_jobs = []

				except subprocess.CalledProcessError, e:
					print "Error processing stdout output: %s for sample %s\n" % (e.output, sample_name)
				except Exception:
					print "Error processing %s. " % sample_name

			idx += 1
	runHtseq(htseq_jobs)


def main():

	global GENOME_DIR
	global GENOME_FASTA
	global SJDB_GTF

	GENOME_DIR="/mnt/genome/STAR/HG38-PLUS/" # change
	GENOME_FASTA="/mnt/genome/hg38-plus/hg38-plus.fa" # change
	SJDB_GTF="/mnt/genome/hg38-plus/hg38-plus.gtf" # change

	ref_genome_file = 'hg38-plus.tgz'
	ref_genome_star_file = 'HG38-PLUS.tgz'

	taxon   = os.environ['TAXON']
	if taxon == 'mus':
		GENOME_DIR="/mnt/genome/STAR/MM10-PLUS/"
		GENOME_FASTA="/mnt/genome/mm10-plus/mm10-plus.fa" # change
		SJDB_GTF="/mnt/genome/mm10-plus/mm10-plus.gtf" # change
		ref_genome_file = 'mm10-plus.tgz'
		ref_genome_star_file = 'MM10-PLUS.tgz'

	# download the genome data

	command = "mkdir -p /mnt/genome; aws s3 cp s3://czi-hca/ref-genome/%s /mnt/genome" % ref_genome_file
	print command
	subprocess.check_output(command, shell=True)

	command = "cd /mnt/genome; tar xvfz %s " % ref_genome_file
	print command
	subprocess.check_output(command, shell=True)

	command = "mkdir -p /mnt/genome/STAR; aws s3 cp s3://czi-hca/ref-genome/STAR/%s /mnt/genome/STAR" % ref_genome_star_file
	print command
	subprocess.check_output(command, shell=True)

	command = "cd /mnt/genome/STAR; tar xvfz %s" % ref_genome_star_file
	print command
	subprocess.check_output(command, shell=True)


	sys.stdout.flush()

	# Load Genome Into Memory
	command = "%s --genomeDir %s --genomeLoad LoadAndExit" % (STAR, GENOME_DIR)
	print command
	subprocess.check_output(command, shell=True)

	# run through the samples
	global S3_BUCKET
	S3_BUCKET = os.environ['S3_BUCKET']

	num_partitions = int(os.environ['NUM_PARTITIONS'])
	partition_id = int(os.environ['PARTITION_ID'])
	doc_ids = os.environ['EXP_IDS'].split(",")
	run(doc_ids, num_partitions, partition_id)

	# Remove Genome Into Memory
	command = "%s --genomeDir %s --genomeLoad Remove" % (STAR, GENOME_DIR)
	print command
	subprocess.check_output(command, shell=True)


if __name__ == "__main__":
	main()
