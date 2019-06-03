#! /usr/bin/env python
import sqlite3
import csv
import os
import sys
import pandas as pd
import random
import re
import tarfile
import argparse


def process(merged_file, output_dir):
    mergedDf = pd.read_csv(merged_file, sep='\t')
    mergedDf = mergedDf.loc[:, ['eGene', 'First_Level']]
    pathway_dict = {}
    pathway_set = set()
    for idx, row in mergedDf.iterrows():
        gene = row['eGene']
        if gene in pathway_dict:
            continue
        gene_dict = {}
        pathways = row['First_Level'].split(',')
        pathways = [pathway.strip() for pathway in pathways]
        for pathway in pathways:
            pathway_set.add(pathway)
            gene_dict[pathway] = 1
        pathway_dict[gene] = gene_dict
    outfile = open(os.path.join(output_dir, 'pathway_set_upset.txt'), 'w')
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(['Gene'] + list(sorted(pathway_set)))
    for gene in pathway_dict:
        gene_pathways = []
        for pathway in sorted(pathway_set):
            try:
                gene_pathways.append(pathway_dict[gene][pathway])
            except KeyError:
                gene_pathways.append(0)
        writer.writerow([gene] + gene_pathways)
    outfile.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--merged_file", required=True,
        help="the merged.txt file created with viz.Rmd")
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Directory in which to write output.")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    #merged_file = '../analysis/merged.txt'
    #output_dir = '../analysis/'
    process(args.merged_file, args.output_dir)
