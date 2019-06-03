#!/usr/bin/env python
import pandas as pd
import numpy as np
import csv
import sys
import os
import numbers
import argparse


"""
Outputs files to used to draw phenograms at PhenoGram from Ritchie Lab
(http://visualization.ritchielab.org/phenograms) of metabolites on karyotypes.

The output files include
    1. phenogram_all.txt, containing all spatial variant-gene pairs.
    2. phenogram_metabolism, only variant-gene pairs involved in metabolism

Files have the ff structure:
    1. SNP
    2. ANNOTATION #Gene name
    3. CHR # Gene chromosome
    4. POS # Gene start
    5. SNP_PHENOTYPE # Biochemical associated with SNP locus in Shin et al 2014
    6. END # Gene end
    7. First_Level # Gene's list of Kegg's super pathway categories
    8. PHENOTYPE # Gene's sub pathway category if First_Level is Metabolism
    9. Metabolism # A boolean field for genes involved in metabolism
"""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e", "--eqtl_file", required=True,
        help="significant_eqtls.txt from CoDeS3D run.")
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Directory to write results.")
    parser.add_argument(
        "-k", "--kegg_map", default="../data/kegg_pathway_map.txt",
        help="A map of KEGG pathways.")
    parser.add_argument(
        "-r", "--kegg_results", required=True,
        help="Results from KEGG pathways run of eGenes in text format.")
    parser.add_argument(
        "-m", "--metab_file", default="../data/shin_2014_suppl_tables.xlsx",
        help="Supplementary spreadsheet from Shin et al.")
    return parser.parse_args()


def parse_kegg_ref(kegg_map):
    print('Parsing Kegg Pathway reference...')
    kegg_ref = {}
    kegg_map_file = open(kegg_map, 'r')
    map_reader = csv.reader(kegg_map_file, delimiter='\t')
    next(map_reader, None)
    first_level = ''
    second_level = ''
    for row in map_reader:
        row = row[0].strip()
        if not row or row[0].strip().startswith('#'):
            continue
        if row[0].strip().startswith('*'):
            first_level = row[1:]
            if not first_level in kegg_ref.keys():
                kegg_ref[first_level] = {}
            continue
        try:
            int(row[0].strip()[0])
            kegg_ref[first_level][second_level].append(row[7:])
        except:
            second_level = row
            if second_level not in kegg_ref[first_level].keys():
                kegg_ref[first_level][second_level] = []
    return kegg_ref


def map_pathways():
    kegg_data = []
    pathways = {}
    kegg_not_found = []
    kegg_f = open(kegg_results, 'r')
    kegg = csv.reader(kegg_f, delimiter='\t')
    next(kegg, None)
    print('Mapping eGenes pathways...')
    for row in kegg:
        if row:
            kegg_data.append(row)
    for item in kegg_data[0][0].split():
        if item.startswith('hsa'):
            kegg_not_found.append(item[4:])
    pathway = ''
    for row in kegg_data[3:]:
        if row[0].endswith(')'):
            try:
                # Last items in pathway sub-headers are bracketed numbers
                int(row[0][-2])
                row1 = row[0].split()
                row1 = row1[1: len(row1) - 5]
                k = ' '.join(row1)
                pathway = k
                pathways[pathway] = []
            except:
                pass
        else:
            row1 = row[0].split()[1]
            pathways[pathway].append(row1[: len(row1) - 1])
    genes = {}
    for p in sorted(pathways.keys()):
        for gene in pathways[p]:
            if not gene in genes.keys():
                genes[gene] = {'pathway': [],
                               'mapping': {}}
            genes[gene]['pathway'].append(p)
    for gene in sorted(genes.keys()):
        for first_level in kegg_ref:
            for second_level in kegg_ref[first_level]:
                if set(genes[gene]['pathway'])\
                   .intersection(set(kegg_ref[first_level][second_level])):
                    if first_level not in genes[gene].keys():
                        genes[gene]['mapping'][first_level] = set([])
                    genes[gene]['mapping'][first_level].add(second_level)
    return genes, pathways


def process():
    print('Processing results...')
    metab_df = pd.read_excel(metab_file, sheet_name="Table S6",
                             skiprows=[0, 1, 2])
    eqtl_df = pd.read_table(eqtl_file, delimiter='\t')
    joined_df = eqtl_df.set_index('SNP').join(metab_df.set_index('SNP'))
    joined_df.rename(columns={'Gene_Chromosome': 'CHR',
                              'Gene_Start': 'POS',
                              'Most associated metabolite or ratio': 'SNP_PHENOTYPE',
                              'Gene_End': 'END',
                              'Gene_Name': 'ANNOTATION'},
                     inplace=True)
    joined_df['rsID'] = joined_df.index
    no_dups = joined_df.loc[:, ['ANNOTATION',
                                'CHR', 'POS', 'SNP_PHENOTYPE', 'END', 'rsID']]
    #print(no_dups.loc[no_dups.duplicated(), :])
    no_dups = no_dups.drop_duplicates(keep='first')
    no_dups['First_Level'] = 'Unknown'
    no_dups['PHENOTYPE'] = ''
    no_dups['Metabolism'] = 'False'
    print(no_dups.columns)
    for idx, annot in enumerate(no_dups['ANNOTATION']):
        # if len(annot) > 10:  # PhenoGram will not accept len(ANNOTATION) > 10
        #    no_dups.iloc[idx, 0] = annot[:8]
        pathway = ''
        gene = no_dups.iloc[idx, 0]
        try:
            first_levels = list(genes[gene]['mapping'].keys())
            if 'Metabolism' in first_levels:
                pathway = 'Metabolism'
                second_level = ', '.join(
                    list(genes[gene]['mapping']['Metabolism']))
                # Make secondary grouping more informative.
                if second_level == 'Global and overview maps':
                    if 'Metabolic pathways' in genes[gene]['pathway']:
                        second_level = 'Other metabolic pathways'
                elif second_level == 'Glycan biosynthesis and metabolism':
                    second_level = 'Carbohydrate metabolism'
                elif second_level == 'Metabolism of other amino acids':
                    second_level = 'Amino acid metabolism'
                no_dups.iloc[idx, 7] = second_level
                no_dups.iloc[idx, 8] = 'True'
            if len(first_levels) > 0:
                pathway = ', '.join(first_levels)
            else:
                pathway = first_levels[0]
            no_dups.iloc[idx, 6] = pathway
        except:
            pass
    met_df = no_dups[no_dups.Metabolism == 'True']
    met_df.to_csv(os.path.join(output_dir, 'phenogram_metabolism.txt'),
                  header=True, sep='\t', mode='w')
    no_dups.to_csv(os.path.join(output_dir, 'phenogram_all.txt'),
                   header=True, sep='\t', mode='w')
    print('Done.')


if __name__ == '__main__':
    args = parse_args()
    eqtl_file = args.eqtl_file
    metab_file = args.metab_file
    kegg_map = args.kegg_map
    kegg_results = args.kegg_results
    output_dir = args.output_dir

    #kegg_file = '../results/kegg_results.txt'
    #kegg_map = '../data/kegg_pathway_map.txt'
    #output_dir = '../results/'
    kegg_ref = parse_kegg_ref(kegg_map)
    genes, pathways = map_pathways()
    process()
