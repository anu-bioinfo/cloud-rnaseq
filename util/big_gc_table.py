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

def taxon_name_to_id(taxon_name):
    if taxon_name == 'mus':
        return MUS_NCBI_TAXON_ID
    if taxon_name == 'homo':
        return HOMO_NCBI_TAXON_ID

def generate_metadata_objects(gse_file):
    gse = GEOparse.get_GEO(filepath=gse_file)
    # output
    # {sample_data: , series_data: , platform_data:}
    series_data = gse.metadata
    platform_data = []
    samples = {}
    for key, obj in gse.gpls.iteritems():
        obj_hash = obj.metadata
        obj_hash['gpl_key'] = key
        platform_data.append(obj_hash)

    for key, obj in gse.gsms.iteritems():
        sample_data = obj.metadata
        samples[key] = {'sample_data':sample_data, 'series_data':series_data,
                        'platform_data':platform_data}
    return samples

def extract_top_metadata_sample_fields(samples_obj, taxon_id, doc_id):
    output = []
    sample_list =[]
    for (key, sample) in samples_obj.iteritems():
        taxons = sample['sample_data']['taxid_ch1']
        if taxon_id not in taxons: # taxon doesn't match
            continue
        sample_list.append(key)
        row = [doc_id]
        for f in METADATA_TOP_SAMPLE_FIELDS:
            val = sample['sample_data'].get(f, '')
            if isinstance(val, list):
                val = ";".join(val)
            row.append(val)
        output.append(row)
    return (sample_list, output)

def extract_top_metadata_series_fields(samples_obj, doc_id):
    output = [doc_id]
    series_data = samples_obj.values()[0]['series_data']
    for f in METADATA_TOP_SERIES_FIELDS:
        val = series_data.get(f, '')
        if isinstance(val, list):
            val = ";".join(val)
        output.append(val)
    return output

def gen_metadata_file(doc_id, output_dir):
    try:
        s3_source = S3_BUCKET + '/' + doc_id + '/metadata/'
        command = "mkdir -p %s;mkdir -p _tmp/%s; aws s3 cp %s _tmp/%s/ --recursive --exclude '*' --include '*family.soft.gz'" % (output_dir, doc_id, s3_source, doc_id)
        print command
        output = subprocess.check_output(command, shell=True)

        command = "find _tmp/%s -maxdepth 1 -name '*family.soft.gz'" % doc_id
        print command
        output = subprocess.check_output(command, shell=True)
        if output:
            meta_files = output.rstrip().split("\n")
            samples_obj = generate_metadata_objects(meta_files[0])
            output_file = "%s/%s.metadata.json" % (output_dir, doc_id)
            with open(output_file, 'wb') as csvfile:
                csvfile.write(json.dumps(samples_obj))
            return samples_obj
    except Exception:
        return


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
                SRR_GSM_MAPPING[row[0]] = row[1]
    except Exception:
        return {}

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
    get_srr_gsm_mapping(doc_id)
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


def generate_gc_table(htseq_list, gene_to_idx_table, doc_sample_list):
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
            cell_info['gsm_id'] = SRR_GSM_MAPPING.get(m.group(2)) or ''

        if len(cell_info['doc_id']) > 0:
            doc_id = cell_info['doc_id']
            sample_list = doc_sample_list[doc_id]
            if len(cell_info['gsm_id']) > 0 and (cell_info['gsm_id'] not in sample_list):
                continue # Skip if sample is not included (due to different taxon)
            doc_to_counts[doc_id] = doc_to_counts.get(doc_id, 0) + 1
        cells.append(cell_info)
    return (doc_to_counts, cells)

def output_logfiles_to_csv(log_list, log_csv, doc_sample_list):
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
            cell_info['gsm_id'] = SRR_GSM_MAPPING.get(m.group(2)) or ''
        if len(cell_info['doc_id']) > 0:
            sample_list = doc_sample_list[cell_info['doc_id']]
            if len(cell_info['gsm_id']) > 0 and (cell_info['gsm_id'] not in sample_list):
                continue # Skip if sample is not included (due to different taxon)
        cells.append(cell_info)
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

def output_metadata_csv(metadata_table, output_file_name, headers):
    headers = ['doc_id'] + headers
    with open(output_file_name, 'wb') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        gc_writer.writerow(headers)
        for row in metadata_table:
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
    parser.add_argument('-f', action="store", dest='output_csv', default=False)
    parser.add_argument('-l', action="store", dest='log_csv', default=False)
    parser.add_argument('-t', action="store", dest='taxon', default=False, help='mus or homo')
    parser.add_argument('-m', action="store", dest='metadata_prefix', default=False,
        help='meta data prefix. output would <metadata_prefix>.tgz')
    results = parser.parse_args()
    if results.s3_path and results.output_csv and results.log_csv and results.taxon and results.metadata_prefix:
        S3_BUCKET=results.s3_path
        taxon_id = taxon_name_to_id(results.taxon)
        htseq_list =[]
        log_list = []
        doc_list = get_doc_list(S3_BUCKET)
        metadata_list = []
        series_summary_list = []
        doc_to_sample_list = {}
        for doc_id in doc_list:
            htseq_list += get_htseq_files_from_s3(doc_id, results.taxon)
            log_list += get_log_files_from_s3(doc_id, results.taxon)
            samples_obj = gen_metadata_file(doc_id, results.metadata_prefix)
            taxon_id = taxon_name_to_id(results.taxon)
            (doc_samples, samples_metadata) = extract_top_metadata_sample_fields(samples_obj,
                    taxon_id, doc_id)
            series_metadata = extract_top_metadata_series_fields(samples_obj, doc_id)
            doc_to_sample_list[doc_id] = doc_samples
            series_summary_list.append(series_metadata)
            metadata_list += samples_metadata

        # generating tar ball for all the metadata
        command = "tar cvfz %s.tgz %s" % (results.metadata_prefix, results.metadata_prefix)
        print command
        output = subprocess.check_output(command, shell=True)

        # creating a temp work space
        command = "mkdir -p _tmp/%s" % doc_id
        print command
        output = subprocess.check_output(command, shell=True)


        # output metadata

        if len(metadata_list) > 0:
            output_file_name = results.metadata_prefix + '.samples.csv'
            output_metadata_csv(metadata_list, output_file_name, METADATA_TOP_SAMPLE_FIELDS)


        if len(htseq_list) > 0:
            # start the mapping
            gene_to_idx_table = get_gene_to_idx_mapping(htseq_list[0])
            (doc_to_counts, gene_cell_table) = generate_gc_table(htseq_list, gene_to_idx_table,
                    doc_to_sample_list)
            output_htseq_to_csv(gene_cell_table, results.output_csv, gene_to_idx_table)
            if len(series_summary_list) > 0:
                output_series_list = []
                for row in series_summary_list:
                    doc_id = row[0]
                    if doc_to_counts.get(doc_id, 0) > 0:
                        output_series_list.append(row)
                output_file_name = results.metadata_prefix + '.series.csv'
                output_metadata_csv(output_series_list, output_file_name, METADATA_TOP_SERIES_FIELDS)

        else:
            print "No htseq files for %s " % (S3_BUCKET)

        if len(log_list) > 0:
            output_logfiles_to_csv(log_list, results.log_csv, doc_to_sample_list)
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

