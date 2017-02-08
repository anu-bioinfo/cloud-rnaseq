import argparse
import subprocess
import os

from Bio import Entrez

import json
import redis
import re

REDIS_STORE = redis.StrictRedis(host='localhost', port=6379, db=2)
FTP_HOST = "ftp://ftp-trace.ncbi.nlm.nih.gov"
DEST_DIR = "/data/hca/"
FASTQ_COMMAND = "/usr/local/bin/sratoolkit/fastq-dump -split-3 "
WGET = '/usr/bin/wget '
SEQUENCE_COMMAND = '/home/yunfang/code/geo_download/scripts/1pass_and_htseq.sh'
MIN_FILE_SIZE = 5000

def run_sample(download_path):
	''' Example: 
	   run_sample("/sra/sra-instant/reads/ByStudy/sra/SRP/SRP055/SRP055440//SRR1813952/SRR1813952.sra" 
	'''
	filename = os.path.basename(download_path)
	sample_name = os.path.splitext(filename)[0]
	dest_file = DEST_DIR + 'rawdata/' + filename
	download_url = FTP_HOST + download_path
	# Download data
	print "Fetch the file through FTP"
	if (not os.path.exists(dest_file)) or (os.path.getsize(dest_file) < MIN_FILE_SIZE):
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
	
	print "Run sequencing"
	command = "ls " + DEST_DIR + 'rawdata/' + sample_name + '*.fastq' 
	fastq_files = subprocess.check_output(command, shell=True).replace("\n", " ")
	result_dir = DEST_DIR + 'results/' + sample_name

	command = SEQUENCE_COMMAND + ' ' + result_dir + ' ' + fastq_files
	print command
	output = subprocess.check_output(command, shell=True)
	print output
	
def run(doc_id):
	''' Example: run("200066217") '''
	
def main():
	parser = argparse.ArgumentParser(
		description='Download data from a GEO dataset, sequence and generate expression counts')


if __name__ == "__main__":
	main()
