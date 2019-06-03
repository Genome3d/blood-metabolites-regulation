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


def process(protein_atlas_fp, eqtl_fp, num_sim, output_dir):
    print('Preparing gene list...')
    global summary
    summary = []
    eqtl_file = pd.read_csv(eqtl_fp, sep='\t')
    genes = eqtl_file['Gene_Name'].unique()
    global protein_atlas
    protein_atlas = pd.read_csv(protein_atlas_fp, compression='zip', sep='\t',
                                keep_default_na=False)
    gene_pool = protein_atlas['Gene'].unique()
    num_genes = len(set(genes).intersection(set(gene_pool)))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    print('Bootstraping...')
    for j in range(num_sim):
        random_gene_list = random.sample(gene_pool, num_genes)
        # for i in random_numbers:
        #    random_gene_list.append(
        #        gene_pool[i])
        simulate(j, random_gene_list, output_dir)
    print('There are {} genes in the Protein Atlas'.format(len(gene_pool)))
    print('{} of your {} genes are annotated in the Protein Atlas'
          .format(num_genes, len(genes)))


def simulate(j, gene_set, output_dir):
    protein_class = {}
    gene_list = protein_atlas[protein_atlas.Gene.isin(gene_set)]
    gene_list = gene_list.drop_duplicates(['Gene'])
    i = len(gene_list)
    for idx, row in gene_list[gene_list.duplicated(['Gene'])].iterrows():
        print(row['Gene'])
    for idx, row in gene_list.iterrows():
        for class_ in row['Protein class'].split(','):
            class_ = class_.strip()
            if not class_ in protein_class:
                protein_class[class_] = []
            protein_class[class_].append(row['Gene'])
    simfile = open(os.path.join(output_dir, str(j) + '.txt'), 'w')
    simwriter = csv.writer(simfile, delimiter='\t')
    for class_ in sorted(protein_class):
        simwriter.writerow((class_, len(protein_class[class_])))
    print('\tset {}:\t{} {} genes'.format(j+1, i, len(gene_set)))
    simfile.close()


def get_protein_class(gene_set):
    """Returns the protein classes of genes

    Args:
    A dictionary of genes:
        {'Gene1': {'name': 'Full gene name', 'comp_ids': [1234, 1235]}
         'Gene2': {'name': 'Full gene name', 'comp_ids': [2345, 2345]}

    Output:
    A dictionary of protein class of gene products:
        {'Gene1': ['Protein family name'],
         'Gene2': ['Protein family name']}
    """
    chemblDB = sqlite3.connect(chembl_file)
    chemblDB.text_factory = str
    cur = chemblDB.cursor()
    protein_class = {}
    for gene in gene_set:
        protein_class[gene] = set()
        for compid in gene_set[gene]['comp_ids']:
            cur.execute(
                """SELECT pref_name FROM protein_classification WHERE protein_class_id =
                (SELECT protein_class_id FROM component_class WHERE component_id = ?)
                """, (compid,))
            """
            aspect values:
            F=molecular function; P=biological process; C=cellular component
            """
            protein_class[gene].add(cur.fetchone()[0])
    cur.close()
    chemblDB.close()
    for gene in protein_class:
        protein_class[gene] = list(protein_class[gene])
    return protein_class


def collate(num_sim, output_dir):
    class_dict = {}
    num_sim = 0
    for sim_file in os.listdir(output_dir):
        aggre_dict = {'Blood group antigen proteins': 0,
                      'Cancer-related genes': 0,
                      'Candidate cardiovascular disease genes': 0,
                      'CD markers': 0,
                      'Citric acid cycle related proteins': 0,
                      'Disease related genes': 0,
                      'Enzymes': 0,
                      'FDA approved drug targets': 0,
                      'G-protein coupled receptors': 0,
                      'Nuclear receptors': 0,
                      'Plasma proteins': 0,
                      'Potential drug targets': 0,
                      'Predicted intracellular proteins': 0,
                      'Predicted membrane proteins': 0,
                      'Predicted secreted proteins': 0,
                      'RAS pathway related proteins': 0,
                      'Ribosomal proteins': 0,
                      'Transcription factors': 0,
                      'Transporters': 0,
                      'Voltage-gated ion channels': 0,
                      'RNA polymerase related proteins': 0}
        num_sim += 1
        sim_file = open(os.path.join(output_dir, sim_file), 'r')
        sim = csv.reader(sim_file, delimiter='\t')
        for row in sim:
            if not row[0] in class_dict:
                class_dict[row[0]] = 0
            class_dict[row[0]] += int(row[1])
            aggre_dict[row[0]] = row[1]
        sim_file.close()
        to_file = [v for k, v in aggre_dict.iteritems()]
    class_file = open(os.path.join(
        '../analysis/', 'bootstrap_protein_atlas_class.txt'), 'w')
    sim = csv.writer(class_file, delimiter='\t')
    sim.writerow(['Class', 'Freq', 'Overlap'])
    for c in sorted(class_dict):
        sim.writerow([c, class_dict[c]/float(num_sim), 'Bootstrap'])
    class_file.close()


def aggregate(num_sim, output_dir):
    to_file = []
    num_sim = 0
    for sim_file in os.listdir(output_dir):
        aggre_dict = {'Blood group antigen proteins': 0,
                      'Cancer-related genes': 0,
                      'Candidate cardiovascular disease genes': 0,
                      'CD markers': 0,
                      'Citric acid cycle related proteins': 0,
                      'Disease related genes': 0,
                      'Enzymes': 0,
                      'FDA approved drug targets': 0,
                      'G-protein coupled receptors': 0,
                      'Nuclear receptors': 0,
                      'Plasma proteins': 0,
                      'Potential drug targets': 0,
                      'Predicted intracellular proteins': 0,
                      'Predicted membrane proteins': 0,
                      'Predicted secreted proteins': 0,
                      'RAS pathway related proteins': 0,
                      'Ribosomal proteins': 0,
                      'Transcription factors': 0,
                      'Transporters': 0,
                      'Voltage-gated ion channels': 0,
                      'RNA polymerase related proteins': 0}
        num_sim += 1
        sim_file = open(os.path.join(output_dir, sim_file), 'r')
        sim = csv.reader(sim_file, delimiter='\t')
        for row in sim:
            # if not row[0] in class_dict:
            #    class_dict[row[0]] = 0
            # class_dict[row[0]] += int(row[1])
            aggre_dict[row[0]] = row[1]
        sim_file.close()
        to_file.append([aggre_dict[k] for k in sorted(aggre_dict)])
    class_file = open(os.path.join(
        '../analysis/', 'bootstrap_protein_atlas_class_aggregate.txt'), 'w')
    sim = csv.writer(class_file, delimiter='\t')
    sim.writerow([k for k in sorted(aggre_dict)])
    sim.writerows(to_file)
    class_file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e", "--eqtl_fp", required=True,
        help="significant_eqtls.txt file from CoDeS3D run")
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Directory to write results.")
    parser.add_argument(
        "-b", "--num_sim", default=1000,
        help="Bootstrap value")
    parser.add_argument(
        "-p", "--protein_atlas_fp", default='../data/proteinatlas.tsv.zip',
        help="Protein Atlas database")
    parser.add_argument(
        "-g", "--gene_ref_fp", default='../../codes3d-v1/lib/gene_reference.bed',
        help="Gene reference bed file")
    return(parser.parse_args())


if __name__ == '__main__':
    args = parse_args()
    #protein_atlas_fp = '../data/proteinatlas.tsv.zip'
    #gene_file = '../../codes3d-v1/lib/gene_reference.bed'
    #eqtl_fp = '../results/codes3d_output/significant_eqtls.txt'
    #output_dir = '../analysis/bootstrap_protein_class/'
    num_sim = int(args.num_sim)
    process(args.protein_atlas_fp, args.eqtl_fp, num_sim, args.output_dir)
    collate(num_sim, args.output_dir)
    aggregate(num_sim, args.output_dir)
