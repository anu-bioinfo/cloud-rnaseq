import argparse
import subprocess
import os
import thread
import time

from Bio import Entrez

import json
import redis
import re

REDIS_STORE = redis.StrictRedis(host='localhost', port=6379, db=2)
FTP_HOST = "ftp://ftp-trace.ncbi.nlm.nih.gov"
DEST_DIR = "/data/hca/"
FASTQ_COMMAND = "/usr/local/bin/sratoolkit/fastq-dump -split-3 "
WGET = '/usr/bin/wget '
SEQUENCE_COMMAND = '/home/yunfang/code/hca/geo_download/scripts/1pass_and_htseq.sh'
MIN_FILE_SIZE = 5000
MAX_THREADS = 3
S3_BUCKET = 's3://czsi-sra'

def run_sample(download_path, s3_prefix):
	''' Example: 
	   run_sample("/sra/sra-instant/reads/ByStudy/sra/SRP/SRP055/SRP055440//SRR1813952/SRR1813952.sra" 
	'''
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
	command = "gzip %s" % fastq_files
	print command
	output = subprocess.check_output(command, shell=True)
	print output

	# move data over to s3
        s3_dest = S3_BUCKET + '/' + s3_prefix + '/rawdata/' + sample_name + '/'
        sample_files =  DEST_DIR + 'rawdata/'+ sample_name + '*.*' 
	command = "s3cmd put %s %s" % (sample_files, s3_dest)

	print command
	output = subprocess.check_output(command, shell=True)
	print output

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

def run_samples(sample_list, doc_id, tn):
	for download_path in sample_list:
		while 1:
			try:
				run_sample(download_path, doc_id)
				break
			except CalledProcessError:
				time.sleep(30)
		#print "%d run_sample(%s)" % (tn, download_path)
		#time.sleep(1)
	
def run(doc_id):
	''' Example: run("200066217") '''
	file_list = json.loads(REDIS_STORE.get('sra:' + doc_id) or "[]")
	thread_queues = []
	for i in range(MAX_THREADS):
		thread_queues.append([])
	t_idx = 0
	for item in file_list:
		t_idx = 0 if t_idx >= MAX_THREADS else t_idx
		thread_queues[t_idx].append(item[1])
		t_idx += 1
		#print "gp.run_sample('%s')" % item[1]
	try:
		idx = 1
		for queue in thread_queues:
			thread.start_new_thread(run_samples, (queue, doc_id, idx))
			time.sleep(30)
			idx += 1
	except:	
		print "Error: unable to start thread"
	
	return thread_queues
	
	
def main():
	parser = argparse.ArgumentParser(
		description='Download data from a GEO dataset, sequence and generate expression counts')


if __name__ == "__main__":
	main()
