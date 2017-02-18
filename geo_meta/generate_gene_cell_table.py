#!/usr/bin/python
import json
import time
import os
import sys

from Bio import Entrez
import redis
import redis
import hashlib
import re
from ftplib import FTP
import codecs
import argparse

import GEOparse
import csv
import subprocess

S3_BUCKET = 's3://czi-hca/data'

def get_srr_gsm_mapping(doc_id):
    ''' This doesn't work. Fix '''
    # Downalod the data
    output = {}
    s3_dir = S3_BUCKET + '/' + doc_id + '/metadata/'
    command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*.soft.gz'" % (doc_id, s3_dir, doc_id)
    print command
    x = subprocess.check_output(command, shell=True)
    accession_file = subprocess.check_output("ls _tmp/%s/*%s*.soft.gz" % (doc_id, doc_id[-5:]), shell=True).rstrip()  
    gse = GEOparse.get_GEO(filepath=accession_file)
    for gsm_id, gsm_obj in gse.gsms.iteritems():
        sra_url = gsm_obj.metadata['relation'][1]
        m = re.match(".*term=(.*)", sra_url)
        if m: 
            srr = m.group(1)
            output[srr] = gsm_id
    #print output
    return output
    
def get_htseq_files_from_s3(doc_id):
    s3_dir = S3_BUCKET + '/' + doc_id + '/results/'
    command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*.htseq-count.txt'" % (doc_id, s3_dir, doc_id)
    print command
    output = subprocess.check_output(command, shell=True)

    command = "ls _tmp/%s/*.htseq-count.txt" % doc_id
    output = subprocess.check_output(command, shell=True)
    if output:
        htseq_files = output.rstrip().split("\n")
        #print htseq_files
        return htseq_files

def get_gene_to_idx_mapping(ht_seq_path):
    print "getting gene to idx mapping: %s" % ht_seq_path
    gene_to_idx = {}
    idx = 0
    with open(ht_seq_path, "r") as ins:
        for line in ins: 
            fields = line.rstrip().split("\t")
            gene_to_idx[fields[0]] = idx
            idx += 1
    return gene_to_idx
    

def generate_gc_table(htseq_list, srr_to_gsm_map, gene_to_idx_table):
    cells = []
    for htseq_file in htseq_list: 
        cell_data = [0] * len(gene_to_idx_table)
        with open(htseq_file, "r") as ins:
            for line in ins: 
                fields = line.rstrip().split("\t")
                gene = fields[0]
                exp_cnt = fields[1]
                idx = gene_to_idx_table[gene]
                cell_data[idx] = exp_cnt
        m = re.match(".*?([^\/\.]*)\.htseq-count.txt", htseq_file)
        cell_info= {'data':cell_data, 'srr_id': '', 'gsm_id': ''}
        if m:
            cell_info['srr_id'] = m.group(1)
            #try:
            #    cell_info['gsm_id'] = srr_to_gsm_map[m.group(1)]
            #except Exception: 
            #    print "Couldn't find gsm mapping for %s." % cell_info['srr_id']
        cells.append(cell_info)
    return cells
        
        
def main():
    parser = argparse.ArgumentParser(
        description='Generate gene cell table')
    parser.add_argument('-d', action="store", dest='gds_id', default=False)
    parser.add_argument('-f', action="store", dest='output_csv', default=False)
    results = parser.parse_args()
    if results.gds_id and results.output_csv:
        doc_id = results.gds_id
	htseq_list = get_htseq_files_from_s3(doc_id)
	if len(htseq_list) > 0:
            # creating a temp work space
            command = "mkdir -p _tmp/%s" % doc_id
            print command
            output = subprocess.check_output(command, shell=True)
            
            # start the mapping
	    srr_to_gsm_map = get_srr_gsm_mapping(doc_id)
	    gene_to_idx_table = get_gene_to_idx_mapping(htseq_list[0])
	    gene_cell_table = generate_gc_table(htseq_list, srr_to_gsm_map, gene_to_idx_table)
            
            with open(results.output_csv, 'wb') as csvfile:
                gc_writer = csv.writer(csvfile, delimiter=',')

                header_row = ['GENENAME']
		for cell in gene_cell_table:
                    header_row.append(cell['srr_id'])
                gc_writer.writerow(header_row)

                #header_row = ['GENENAME']
		#for cell in gene_cell_table:
                #    header_row.append(cell['gsm_id'])
                #gc_writer.writerow(header_row)
		# write actual data
		for gene in sorted(gene_to_idx_table.keys()):
                    row = [gene]
                    idx = gene_to_idx_table[gene]
                    for cell in gene_cell_table:
                        row.append(cell['data'][idx])
                    gc_writer.writerow(row)
            # Clean up tmp work space
            command = "rm -rf _tmp/%s" % doc_id
            print command
            output = subprocess.check_output(command, shell=True)
        else:
            print "No htseq files for %s in %s" % (doc_id, S3_BUCKET + '/' + doc_id + '/results/')
            sys.exit(1)
        
    else:
        parser.print_help()
        sys.exit(1)
	
if __name__ == "__main__":
    main()

