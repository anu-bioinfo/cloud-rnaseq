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
SRR_GSM_MAPPING = {} # TODO(yf): make this less hacky
MUS_NCBI_TAXON_ID = '10090'
HOMO_NCBI_TAXON_ID = '9606'
METADATA_TOP_SAMPLE_FIELDS = ['geo_accession', 'title', 'treatment_protocol_ch1',
                     'taxid_ch1', 'characteristics_ch1', 'platform_id']
METADATA_TOP_SERIES_FIELDS = ['title', 'overall_design', 'summary', 'platform_id',
                     'pubmed_id', 'sample_taxid', 'relation']

def generate_sample_to_plate_mapping(s3_path, doc_id):
    s3_source = s3_path + '/' + doc_id + '/rawdata/'
    # Generate mapping
    mapping = {}

    # Get sample list
    command = "aws s3 ls %s --recursive | grep fastq.gz" % s3_source
    print command
    output = subprocess.check_output(command, shell=True).split("\n")
    for fname in output:
        m = re.search("/rawdata/([^/]*)/([^/]*)", fname)
        if m:
            sample = m.group(1)
            fastq_name = m.group(2)
            proper_name = fastq_name.split("__")[0]
            fields = proper_name.split("-")
            if len(fields) > 1:
                mapping[sample] = fields[1]
            else:
                mapping[sample] = 'unknown'
    return mapping

def get_doc_list(s3_path):
    command = "aws s3 ls %s/" % s3_path
    print command
    output = subprocess.check_output(command, shell=True)
    doc_list = []
    if output:
        dir_list = output.rstrip().split("\n")
        for f in dir_list:
            m = re.match(".*?([^\/ ]*)/", f)
            if m and len(m.group(1)) > 0:
                doc_list.append(m.group(1))
    return doc_list

def get_htseq_files_from_s3(doc_id, taxon):
    s3_dir = S3_BUCKET + '/' + doc_id + '/results/'
    try:
        command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*.%s.htseq-count.txt'" % (doc_id, s3_dir, doc_id, taxon)
        print command
        output = subprocess.check_output(command, shell=True)
        command = 'find _tmp/%s -maxdepth 1 -name "*.%s.htseq-count.txt"' % (doc_id, taxon)
        output = subprocess.check_output(command, shell=True)
        if output:
            htseq_files = output.rstrip().split("\n")
            #print htseq_files
            return htseq_files
        else:
            return []
    except Exception:
        print "%s doesn't have any results" % doc_id
        return []

def get_log_files_from_s3(doc_id, taxon):
    s3_dir = S3_BUCKET + '/' + doc_id + '/results/'
    try:
        command = "mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*.%s.log.final.out'" % (doc_id, s3_dir, doc_id, taxon)
        print command
        output = subprocess.check_output(command, shell=True)

        command = 'find _tmp/%s -maxdepth 1 -name "*.%s.log.final.out"' % (doc_id, taxon)
        output = subprocess.check_output(command, shell=True)
        if output:
            log_files = output.rstrip().split("\n")
            #print log_file
            return log_files
        else:
            return []
    except Exception:
        print "%s doesn't have any results" % doc_id
        return []

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

def generate_gc_table(htseq_list, gene_to_idx_table, sample_to_plate_map):
    cells = []
    doc_to_counts = {}
    for htseq_file in htseq_list:
        cell_data = [0] * len(gene_to_idx_table)
        with open(htseq_file, "r") as ins:
            for line in ins:
                fields = line.rstrip().split("\t")
                gene = fields[0]
                exp_cnt = fields[1]
                idx = gene_to_idx_table[gene]
                cell_data[idx] = exp_cnt
        m = re.match(".*?/([^/]*)/([^\/ ]*)\.(.*)\.htseq-count.txt", htseq_file)
        cell_info= {'data':cell_data, 'srr_id': '', 'gsm_id': '', 'taxon': '', 'doc_id': ''}
        if m:
            cell_info['doc_id'] = m.group(1)
            cell_info['srr_id'] = m.group(2)
            cell_info['taxon'] = m.group(3)
            cell_info['gsm_id'] = sample_to_plate_map.get(m.group(2)) or ''

        if len(cell_info['doc_id']) > 0:
            doc_id = cell_info['doc_id']
            doc_to_counts[doc_id] = doc_to_counts.get(doc_id, 0) + 1
        cells.append(cell_info)
    return (doc_to_counts, cells)

def generate_log_table(log_list, sample_to_plate_map):
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
        m = re.match(".*?/([^/]*)/([^\/ ]*)\.(.*)\.log.final.out", log_file)
        cell_info= {'data':cell_data, 'srr_id': '', 'gsm_id': '', 'taxon': '', 'doc_id': ''}
        if m:
            cell_info['doc_id'] = m.group(1)
            cell_info['srr_id'] = m.group(2)
            cell_info['taxon'] = m.group(3)
            cell_info['gsm_id'] = sample_to_plate_map.get(m.group(2)) or ''
        cells.append(cell_info)
    return (headers, cells)

def output_logs_to_csv(log_table, log_csv, headers):
    cells = log_table
    with open(log_csv, 'wb') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        header_row = ['FIELD_NAMES']
        for cell in cells:
            header_row.append(cell['srr_id'])
        gc_writer.writerow(header_row)

        header_row = ['EXP_ID']
        for cell in cells:
            header_row.append(cell['doc_id'])
        gc_writer.writerow(header_row)

        header_row = ['TAXON']
        for cell in cells:
            header_row.append(cell['taxon'])
        gc_writer.writerow(header_row)

        header_row = ['GSM']
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

        header_row = ['EXP_ID']
        for cell in gene_cell_table:
            header_row.append(cell['doc_id'])
        gc_writer.writerow(header_row)

        header_row = ['TAXON']
        for cell in gene_cell_table:
            header_row.append(cell['taxon'])
        gc_writer.writerow(header_row)

        header_row = ['GSM']
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
    global S3_BUCKET

    parser = argparse.ArgumentParser(
        description='Generate gene cell table')
    parser.add_argument('-s', action="store", dest='s3_path', default=False)
    parser.add_argument('-d', action="store", dest='output_dir', default=False, help='output dir')
    parser.add_argument('-t', action="store", dest='taxon', default=False, help='taxon')
    results = parser.parse_args()
    if results.s3_path and results.output_dir and results.taxon:

        command = "mkdir -p %s" % results.output_dir
        print command
        output = subprocess.check_output(command, shell=True)

        S3_BUCKET=results.s3_path
        htseq_list =[]
        log_list = []
        doc_list = get_doc_list(S3_BUCKET)
        metadata_list = []
        series_summary_list = []
        doc_to_sample_list = {}
        sample_to_plate_mapping = {}
        for doc_id in doc_list:
            sample_to_plate_mapping.update(generate_sample_to_plate_mapping(results.s3_path, doc_id))
            htseq_list += get_htseq_files_from_s3(doc_id, results.taxon)
            log_list += get_log_files_from_s3(doc_id, results.taxon)

            # creating a temp work space
            command = "mkdir -p _tmp/%s" % doc_id
            print command
            output = subprocess.check_output(command, shell=True)

        if len(htseq_list) > 0:
            # start the mapping
            gene_to_idx_table = get_gene_to_idx_mapping(htseq_list[0])
            (doc_to_counts, gene_cell_table) = generate_gc_table(htseq_list, gene_to_idx_table,
                    sample_to_plate_mapping)
            # break the data based on plate id
            plate_to_cells = {}
            for cell in gene_cell_table:
                plate_id = cell['gsm_id']
                plate_to_cells[plate_id] = plate_to_cells.get(plate_id, []) + [cell]
            for plate_id, cells in plate_to_cells.iteritems():
                output_htseq_to_csv(cells, results.output_dir+'/'+plate_id+'.htseq-count.csv', gene_to_idx_table)
        else:
            print "No htseq files for %s " % (S3_BUCKET)

        if len(log_list) > 0:
            plate_to_cells = {}
            (headers, log_cell_table) = generate_log_table(log_list, sample_to_plate_mapping)
            for cell in log_cell_table:
                plate_id = cell['gsm_id']
                plate_to_cells[plate_id] = plate_to_cells.get(plate_id, []) + [cell]
            for plate_id, cells in plate_to_cells.iteritems():
                output_logs_to_csv(cells, results.output_dir+'/'+plate_id+'.log.csv', headers)
        else:
            print "No log files for %s" % (S3_BUCKET)

        # Clean up tmp work space
        for doc_id in doc_list:
            command = "rm -rf _tmp/%s" % doc_id
            print command
            output = subprocess.check_output(command, shell=True)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

