#!/usr/bin/env python
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
    # Downalod the data
    try:
        output = {}
        s3_source = S3_BUCKET + '/' + doc_id + '/metadata/srr_to_gsm.csv'
        command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ " % (doc_id, s3_source, doc_id)
        print command
        x = subprocess.check_output(command, shell=True)

        local_file = "_tmp/%s/srr_to_gsm.csv" % doc_id
        with open(local_file, 'rb') as csvfile:
            csvr = csv.reader(csvfile, delimiter=',')
            for row in csvr:
                output[row[0]] = row[1]
        return output
    except Exception:
        return {}

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

def get_log_files_from_s3(doc_id):
    s3_dir = S3_BUCKET + '/' + doc_id + '/results/'
    command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*.log.final.out'" % (doc_id, s3_dir, doc_id)
    print command
    output = subprocess.check_output(command, shell=True)

    command = "ls _tmp/%s/*.log.final.out" % doc_id
    output = subprocess.check_output(command, shell=True)
    if output:
        log_files = output.rstrip().split("\n")
        #print htseq_files
        return log_files

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
        m = re.match(".*?([^\/ ]*)\.htseq-count.txt", htseq_file)
        cell_info= {'data':cell_data, 'srr_id': '', 'gsm_id': ''}
        if m:
            cell_info['srr_id'] = m.group(1)
            cell_info['gsm_id'] = srr_to_gsm_map.get(m.group(1)) or ''
        cells.append(cell_info)
    return cells

def output_logfiles_to_csv(log_list, log_csv, srr_gsm_map):
    cells = []
    headers = []
    first_file = True
    for log_file in log_list:
        cell_data =[]
        with open(log_file, "r") as ins:
            for line in ins:
                fields = line.rstrip().split("|")
                if first_file:
                    headers.append(fields[0])

                if len(fields) > 1:
                    cell_data.append(fields[1])
                else:
                    cell_data.append("")
        first_file = False
        m = re.match(".*?([^\/ ]*)\.log.final.out", log_file)
        cell_info= {'data':cell_data, 'srr_id': '', 'gsm_id': ''}
        if m:
            cell_info['srr_id'] = m.group(1)
            cell_info['gsm_id'] = srr_gsm_map.get(m.group(1)) or ''
        cells.append(cell_info)
    with open(log_csv, 'wb') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        header_row = ['FIELD_NAMES']
	for cell in cells:
            header_row.append(cell['srr_id'])
        gc_writer.writerow(header_row)

        header_row = ['GENENAME']
	for cell in cells:
            header_row.append(cell['gsm_id'])
        gc_writer.writerow(header_row)
	# write actual data
        for idx in range(len(headers)):
            row = [headers[idx]]
            for cell in cells:
                row.append(cell['data'][idx])
            gc_writer.writerow(row)

def output_htseq_to_csv(gene_cell_table, csv_f, gene_to_idx_table):
    with open(csv_f, 'wb') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        header_row = ['GENENAME']
	for cell in gene_cell_table:
            header_row.append(cell['srr_id'])
        gc_writer.writerow(header_row)

        header_row = ['GENENAME']
	for cell in gene_cell_table:
            header_row.append(cell['gsm_id'])
        gc_writer.writerow(header_row)
	# write actual data
        for gene in sorted(gene_to_idx_table.keys()):
            row = [gene]
            idx = gene_to_idx_table[gene]
            for cell in gene_cell_table:
                row.append(cell['data'][idx])
            gc_writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description='Generate gene cell table')
    parser.add_argument('-s', action="store", dest='s3_path', default=False)
    parser.add_argument('-d', action="store", dest='gds_id', default=False)
    parser.add_argument('-f', action="store", dest='output_csv', default=False)
    parser.add_argument('-l', action="store", dest='log_csv', default=False)
    results = parser.parse_args()
    if results.gds_id and results.output_csv and results.log_csv:
        doc_id = results.gds_id
	htseq_list = get_htseq_files_from_s3(doc_id)
        log_list = get_log_files_from_s3(doc_id)
	global S3_BUCKET
	if results.s3_path:
		S3_BUCKET=results.s3_path

        # creating a temp work space
        command = "mkdir -p _tmp/%s" % doc_id
        print command
        output = subprocess.check_output(command, shell=True)

	if len(htseq_list) > 0:

            # start the mapping
	    srr_to_gsm_map = get_srr_gsm_mapping(doc_id)
	    gene_to_idx_table = get_gene_to_idx_mapping(htseq_list[0])
	    gene_cell_table = generate_gc_table(htseq_list, srr_to_gsm_map, gene_to_idx_table)
            output_htseq_to_csv(gene_cell_table, results.output_csv, gene_to_idx_table)

        else:
            print "No htseq files for %s in %s" % (doc_id, S3_BUCKET + '/' + doc_id + '/results/')

        if len(log_list) > 0:
            output_logfiles_to_csv(log_list, results.log_csv, srr_to_gsm_map)
        else:
            print "No log files for %s in %s" % (doc_id, S3_BUCKET + '/' + doc_id + '/results/')

        # Clean up tmp work space
        command = "rm -rf _tmp/%s" % doc_id
        print command
        output = subprocess.check_output(command, shell=True)

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

