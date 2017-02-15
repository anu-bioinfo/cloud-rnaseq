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
FASTQ_COMMAND = "fastq-dump -split-3 "
WGET = '/usr/bin/wget '
MIN_FILE_SIZE = 5000
MAX_THREADS = 3
S3_BUCKET = 's3://czi-hca/data'
MAX_RETRIES = 3

def run_sample(download_path, s3_prefix):
	''' Example:
	   run_sample("/sra/sra-instant/reads/ByStudy/sra/SRP/SRP055/SRP055440//SRR1813952/SRR1813952.sra"
	'''
	subprocess.check_output("mkdir -p %s" % DEST_DIR + 'rawdata/', shell=True)
	filename = os.path.basename(download_path)
	sample_name = os.path.splitext(filename)[0]
	dest_file = DEST_DIR + 'rawdata/' + filename
	done_file = DEST_DIR + 'rawdata/' + sample_name + '.done'
	download_url = FTP_HOST + download_path
	# Download data
	print "Run sample %s" % sample_name
	print "Fetch the file through FTP"
	if (os.path.exists(done_file)):
		return "%s already downloaded" % download_path

	command = "%s -O %s %s" % (WGET, dest_file, download_url)
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	# sra => fastq
	print "Run fastq-dump"
	command = "cd %s; %s %s" % (DEST_DIR + 'rawdata/', FASTQ_COMMAND, dest_file)
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	# run sequencing and htseq
	command = "ls " + DEST_DIR + 'rawdata/' + sample_name + '*.fastq'
	fastq_files = subprocess.check_output(command, shell=True).replace("\n", " ")

	#compressed the files
	command = "gzip -f %s" % fastq_files
	print command
	output = subprocess.check_output(command, shell=True)
	print output

	# move data over to s3
        sample_files =  DEST_DIR + 'rawdata/'+ sample_name + '*.*'
        file_list = subprocess.check_output("ls %s" % sample_files, shell=True).split("\n")
        s3_dest = S3_BUCKET + '/' + s3_prefix + '/rawdata/' + sample_name + '/'
	for f in file_list:
		if len(f) > 5:
			command = "aws s3 cp %s  %s" % (f, s3_dest)
			print command
			output = subprocess.check_output(command, shell=True)

	# cleanup and mark
	command = "rm -rf %s; if [ $? -eq 0 ]; then touch %s; fi" % (sample_files, done_file)
	print command
	output = subprocess.check_output(command, shell=True)
	print output

	'''
	print "Run sequencing"
	result_dir = DEST_DIR + 'results/' + sample_name

	command = SEQUENCE_COMMAND + ' ' + result_dir + ' ' + fastq_files
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	'''

def run(doc_id, num_partitions, partition_id):
	print "Running partition %d of %d for doc %s" % (partition_id, num_partitions, doc_id)
	file_list = json.loads(REDIS_STORE.get('sra:' + doc_id) or "[]")
	idx = 0
	for sra_item in file_list:
		if idx % num_partitions == partition_id:
			tries = 1
			while True:
				try:
					sra_download_path = sra_item[1]
					run_sample(sra_download_path, doc_id)
					print "%s : %s " % (doc_id, sra_download_path)
					time.sleep(5)
					break
				except subprocess.CalledProcessError as e:
					if tries >= MAX_RETRIES:
						print sra_download_path + " doesn't exist."
						break;
					else:
						print "Error downloading %s. Retry in 5 seconds"
						time.sleep(30*tries)
						tries += 1
				
		idx += 1


def main():
	# copy the redis database and start the redis-server
	rdb_s3_path = os.environ['RDB_S3_PATH'] 
	command = "mkdir -p /mnt/redis; aws s3 cp %s /mnt/redis/" % rdb_s3_path
	print command
	throw_away = subprocess.check_output(command, shell = True)

	subprocess.Popen("cd /mnt/redis/; gunzip -f dump.rdb.gz; redis-server --port %d" % REDIS_PORT, shell=True)
	time.sleep(10)

	doc_id = os.environ['GDS_ID']
	num_partitions = int(os.environ['NUM_PARTITIONS'])
	partition_id = int(os.environ['PARTITION_ID'])
	run(doc_id, num_partitions, partition_id)




if __name__ == "__main__":
	main()
