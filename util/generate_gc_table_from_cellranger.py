#!/usr/bin/env python
import json
import time
import os
import sys

import redis
import hashlib
import re
import codecs
import argparse

import csv
import subprocess

import numpy as np

GENE_MAPPING_FILE = 'genes.tsv'
BARCODE_MAPPING_FILE = 'barcodes.tsv'
GC_MATRIX_FILE = 'matrix.mtx'

def generate_cell_to_genes_map(cellranger_dir):
    matrix_file = cellranger_dir + '/' + GC_MATRIX_FILE
    output = {}
    total_gene_count = 0
    with open(matrix_file, 'rb') as mf:
        #skip first two lines
        mf.readline()
        mf.readline()
        (total_gene_count, total_barcodes, total_molecules) = mf.readline().rstrip().split(' ')
        total_gene_count = int(total_gene_count)
        for line in mf:
            (gene_id, barcode_id, count) = line.rstrip().split(' ')
            barcode_id = int(barcode_id)
            gene_id = int(gene_id)
            count = int(count)
            output[barcode_id] = output.get(barcode_id, []) + [(gene_id, count)]
    return output

def load_gene_name_array(cellranger_dir):
    gene_map_file = cellranger_dir + '/' + GENE_MAPPING_FILE
    output = []
    with open(gene_map_file, 'rb') as gmf:
        for line in gmf:
            (gene_code, gene_name) = line.rstrip().split("\t")
            output.append(gene_name)
    return output

def output_gene_cell_table(output_csv, min_gene_count, cell_to_genes_map, gene_name_array):
    with open(output_csv, 'wb') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        gc_writer.writerow(['cell_id'] + gene_name_array) # header row is the gene name array
        for (cell_id, gene_array) in cell_to_genes_map.iteritems():
            if len(gene_array) >= min_gene_count: # check if the cell_id has sufficient data
                row = list(np.zeros(len(gene_name_array)))
                for count_tuple in gene_array:
                    (gene_id, count) = count_tuple
                    row[gene_id - 1] = count # index from 1
                gc_writer.writerow([cell_id] + row)

def main():
    parser = argparse.ArgumentParser(
        description='Generate gene cell table from cellranger ')
    parser.add_argument('-d', action="store", dest='cellranger_dir', default=False,
            help = 'cellranger output directory that include matrix.mtx and etc.')
    parser.add_argument('-f', action="store", dest='output_csv', default=False,
            help = 'output location for the cell gene table')
    parser.add_argument('-m', action="store", dest='min_gene_count', default=1,
            help = 'minimal unique gene count to include a cell barcode')
    results = parser.parse_args()

    if results.cellranger_dir and results.output_csv:
        gene_name_array = load_gene_name_array(results.cellranger_dir)
        cell_to_genes_map = generate_cell_to_genes_map(results.cellranger_dir)
        output_gene_cell_table(results.output_csv, int(results.min_gene_count), cell_to_genes_map, gene_name_array)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()

