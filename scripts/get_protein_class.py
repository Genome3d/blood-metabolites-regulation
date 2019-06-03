#! /usr/bin/env python
import sqlite3
import csv
import os
import sys
import pandas as pd
import argparse
import configparser


def process(eqtl_file, protein_atlas_fp, output_dir):
    print('Parsing input file...')
    eqtl_file = open(eqtl_file, 'r')
    eqtl_df = pd.read_csv(eqtl_file, sep='\t')
    gene_set = get_gene_desc(eqtl_df['Gene_Name'].unique(), chembl_file)
    get_protein_atlas_class(eqtl_df['Gene_Name'].unique(), protein_atlas_fp)


def get_protein_atlas_class(gene_set, protein_atlas_fp):
    protein_atlas = pd.read_csv(protein_atlas_fp, compression='zip', sep='\t',
                                keep_default_na=False)
    gene_list = []
    protein_class = {}
    i = 0
    for idx, row in protein_atlas.iterrows():
        if row['Gene'] in gene_set and not row['Gene'] in gene_list:
            i += 1
            # gene_list.append(row['Gene'])
            row['Subcellular location'] = row['Subcellular location'].replace(
                '<br>', ', ')
            gene_list.append(
                [row['Gene'],
                 row['Protein class'],
                 row['Subcellular location'],
                 row['RNA tissue category'],
                 row['RNA TS'],
                 row['RNA TS TPM'],
                 row['RNA cell line category'],
                 row['RNA CS'],
                 row['RNA CS TPM']]
            )
            for class_ in row['Protein class'].split(','):
                class_ = class_.strip()
                if not class_ in protein_class:
                    protein_class[class_] = []
                protein_class[class_].append(row['Gene'])

    class_file = open('../analysis/protein_atlas_class.txt', 'w')
    writer = csv.writer(class_file, delimiter='\t')
    writer.writerow(['Class', 'Freq'])
    for class_ in sorted(protein_class):
        writer.writerow((class_, len(protein_class[class_])))
    print(i)
    class_file.close()

    protein_file = open('../analysis/protein_atlas_genes.txt', 'w')
    writer = csv.writer(protein_file, delimiter='\t')
    writer.writerow(['Gene', 'Protein class', 'Subcellular location',
                     'RNA tissue category', 'RNA TS', 'RNA TS TPM',
                     'RNA cell line category', 'RNA CS', 'RNA CS TPM'])
    writer.writerows(gene_list)
    protein_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True,
                        help="Significant_eqtls.txt file from CODeS3D. Or a txt file containing a " +
                        "list of genes with \'Gene_Name\' as the column header.")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/config.conf"),
                        help="The configuration file to be used for this " +
                        "ChEMBL run (default: config.conf)")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="The directory in which to output results")
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    data_dir = os.path.join(os.path.dirname(__file__),
                            config.get("Defaults", "DATA_DIR"))
    protein_atlas_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "PROTEIN_ATLAS_FP"))
    process(args.input, protein_atlas_fp, args.output_dir)
