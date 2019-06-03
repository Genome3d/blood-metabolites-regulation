#!/bin/env python

import pandas as pd
import argparse
from Bio import Entrez, Medline
import time
import csv
import sys
import os
import re
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


def get_gene_synonyms(gene_file, gene_synonyms_file):
    synonymsDf = pd.read_csv(gene_synonyms_file, sep='\t')
    geneDf = pd.read_csv(gene_file, sep='\t')
    geneDf = geneDf[['eGene']]
    geneDf.columns = ['Gene']
    mergedDf = pd.merge(geneDf, synonymsDf, how='left', on='Gene')
    mergedDf = mergedDf.fillna('')
    # mergedDf.to_csv('../analysis/gene_synonyms.txt', sep='\t', index=False)
    mergedDf['Synonymns'] = mergedDf['Synonymns']\
        .apply(lambda x: ', '.join(x.split('|')))
    mergedDf['Synonymns'] = mergedDf['Synonymns'].apply(
        lambda x: x.replace(' ', '" "'))
    mergedDf = mergedDf.sort_values(by=['Gene'])
    # Include synonyms and full gene names
    # mergedDf['Query'] = mergedDf[['Gene', 'Synonymns', 'Description']]\
    #      .apply(lambda x: '{}, {}, {}'.format(x[0], x[1], x[2]), axis=1)
    mergedDf['Query'] = mergedDf[['Gene', 'Synonymns', 'Description']]\
        .apply(lambda x: '"{}, "{}", "{}",'.format(x[0], x[1], x[2]), axis=1)
    # Exclude full gene names
    # mergedDf['Query'] = mergedDf[['Gene', 'Synonymns']]\
    #      .apply(lambda x: '"{}, "{}",'.format(x[0], x[1]), axis=1)
    # Remove dashes from genes without synonymns
    mergedDf['Query'] = mergedDf['Query'].apply(
        lambda x: x.replace(' "-",', ''))
    mergedDf['Query'] = mergedDf['Query'].apply(
        lambda x: x.replace(',', '"[Title/Abstract] OR '))
    mergedDf['Query'] = mergedDf['Query'].apply(
        lambda x: x.replace('OR "', 'OR'))
    mergedDf['Query'] = mergedDf['Query'].apply(
        lambda x: x.replace('""', '"'))
    all_search_results = []

    print('Searching Entrez for genes...')
    for index, gene in mergedDf.iterrows():
        print('\t{}'.format(gene['Gene']))
        all_search_results.append({
            'term': gene['Gene'],
            'from_search': search_entrez(gene['Query'])})

    return(mergedDf, all_search_results)


def parse_phenogram(phenogram_file, all_search_results):
    phenogramDf = pd.read_csv(phenogram_file, sep='\t')
    phenogramDf = phenogramDf.sort_values(by=['eGene'])
    phenogramDf['snp_query'] = phenogramDf['SNP']\
        .apply(lambda x: '{}[Title/Abstract]'.format(x))
    print('Searching Entrez for the SNPs...')
    for index, row in phenogramDf[['SNP', 'snp_query']].drop_duplicates().iterrows():
        print('\t{}'.format(row['SNP']))
        all_search_results.append({
            'term': row['SNP'],
            'from_search': search_entrez(row['snp_query'])})

    phenogramDf['PHENOTYPE'] = phenogramDf['SNP_PHENOTYPE']\
        .apply(lambda x: str(x).replace('--', '/'))
    phenogramDf['PHENOTYPE'] = phenogramDf['PHENOTYPE']\
        .apply(lambda x: str(x).replace('*', ''))
    phenogramDf['metabolite_query'] = phenogramDf['PHENOTYPE']\
        .apply(lambda x: str(x).split('/'))
    # phenogramDf['metabolite_query'] = phenogramDf['metabolite_query']\
    #    .apply(lambda x: ' OR '.join(['{}{}{}{}'.format(
    #        '"', j.strip(), '"', '[Title/Abstract]') for j in x]))
    phenogramDf['metabolite_query'] = phenogramDf['metabolite_query']\
        .apply(lambda x: ' OR '.join(['{}{}{}{}'.format(
            '"', j.strip(), '"', '[Title/Abstract]') for j in x]))

    print('Searching Entrez for the metabolites...')
    for index, row in phenogramDf[['PHENOTYPE', 'metabolite_query']]\
            .drop_duplicates().iterrows():
        print('\t{}'.format(row['PHENOTYPE']))
        all_search_results.append({
            'term': row['PHENOTYPE'],
            'from_search': search_entrez(row['metabolite_query'])})

    return(phenogramDf, all_search_results)


def search_gene_metabolite(mergedDf, phenogramDf, all_search_results, outfile):
    print('Preparing queries...')
    for index, row in mergedDf.iterrows():
        try:
            i = phenogramDf.loc[phenogramDf['eGene'] == row['Gene']].index[0]
        except IndexError:
            print(row['Gene'])
        query = '({}) AND ({})'.format(
            row['Query'],
            phenogramDf['metabolite_query'][i])
        # query = query.replace('Title/Abstract', 'All Fields')
        all_search_results.append({
            'term': (row['Gene'], phenogramDf['PHENOTYPE'][i]),
            'from_search': search_entrez(query)})
    outfile = open(outfile, 'w')
    outwriter = csv.writer(outfile, delimiter='\t')
    outwriter.writerow(
        ['Gene', 'Metabolite',  'Common_Ids', 'Common_Ids_List'])
    print("Processing Entrez results...")
    for result in all_search_results:
        id_list = []
        for row in result['from_search']:
            id_list += row['IdList']
        outwriter.writerow([result['term'][0],
                            result['term'][1],
                            len(id_list),
                            list(set(id_list))])
    outfile.close()

    return all_search_results


def search_snp_metabolite(mergedDf, phenogramDf, all_search_results, outfile):
    print('Preparing queries...')
    for index, row in mergedDf.iterrows():
        try:
            i = phenogramDf.loc[phenogramDf['eGene'] == row['Gene']].index[0]
        except IndexError:
            print(row['Gene'])
        query = '({}) AND ({})'.format(
            row['Query'],
            phenogramDf['metabolite_query'][i])
        # query = query.replace('Title/Abstract', 'All Fields')
        all_search_results.append({
            'term': (row['Gene'], phenogramDf['PHENOTYPE'][i]),
            'from_search': search_entrez(query)})

    outfile = open(outfile, 'w')
    outwriter = csv.writer(outfile, delimiter='\t')
    outwriter.writerow(
        ['Gene', 'Metabolite',  'Common_Ids', 'Common_Ids_List'])
    print("Processing Entrez results...")
    for result in all_search_results:
        id_list = []
        for row in result['from_search']:
            id_list += row['IdList']
        outwriter.writerow([result['term'][0],
                            result['term'][1],
                            len(id_list),
                            list(set(id_list))])
    outfile.close()

    return all_search_results


def write_search_results(search_results, outfile):
    outfile = open(outfile, 'w')
    outwriter = csv.writer(outfile, delimiter='\t')
    processed_results = {}
    print("Processing Entrez results...")
    for result in search_results:
        id_list = []
        for row in result['from_search']:
            # print(row.keys())
            id_list += row['IdList']
        processed_results[result['term']] = list(set(id_list))
    outwriter.writerow(
        ['Term1', 'Term2', 'Term1_Ids', 'Term2_Ids', 'Common_Ids', 'Common_Ids_List'])
    for term1 in processed_results:
        for term2 in processed_results:
            common_ids = list(set(processed_results[term1])
                              & set(processed_results[term2]))

            outwriter.writerow(
                (term1, term2,
                 len(processed_results[term1]),
                 len(processed_results[term2]),
                 len(common_ids),
                 ', '.join(common_ids)))
    outfile.close()


def query_entrez():
    info_handle = Entrez.einfo(db="pubmed")
    entrez_info = Entrez.read(info_handle)
    info_handle.close()
    '''
    TITL Title
    WORD Text Word
    TIAB Title/Abstract
    '''
    to_search = ['FADS1', 'fatty acid desaturase', 'rs']
    # search_term = "TDO2[Title/Abstract] OR tryptophan 2,3-dioxygenase[Title/Abstract] OR HYPTRP[Title/Abstract]"
    search_term = "Indoleamine 2,3-dioxygenase 1[Title/Abstract]"
    search_entrez(search_term)
    '''
    for field in entrez_info['DbInfo']['FieldList']:
        print('%(Name)s\t %(FullName)s\t %(Description)s' % field)
        # print(field['Name'], field['Description'])
    '''


def search_entrez(search_term):
    search_results = []
    max_return = 100000
    attempt = 1
    while attempt <= 3:
        try:
            handle = Entrez.esearch(db="pubmed", term=search_term, usehistory="y",
                                    retmax=100000, sort="pub+date")
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server {}".format(err.code))
                print("Attempt {} of 3".format(attempt))
                attempt += 1
                time.sleep(15)
            else:
                raise
    result = Entrez.read(handle)
    handle.close()
    count = int(result['Count'])
    # print('\t\t{}'.format( count))
    search_results.append(result)
    batches = 1
    if result['Count'] > max_return:
        batches = count / max_return
        if count % max_return > 0:
            batches += 1
    if count > max_return:
        for batch in range(1, batches):
            attempt = 1
            while attempt <= 3:
                # print(batch)
                try:
                    handle = Entrez.esearch(db="pubmed",
                                            term=search_term,
                                            usehistory="y", retmax=max_return,
                                            retstart=max_return * batch)
                    break
                except HTTPError as err:
                    if 500 <= err.code <= 599:
                        print("Received error from server {}".format(err.code))
                        print("Attempt {} of 3".format(attempt))
                        attempt += 1
                        time.sleep(15)
                    else:
                        raise
            search_results.append(Entrez.read(handle))
            handle.close()
    return search_results


def download_articles(assoc_file_fp):
    assoc_file = open(assoc_file_fp, 'r')
    reader = csv.reader(assoc_file, delimiter='\t')
    next(reader, None)
    i = 0
    p = re.compile(r'\d+')
    print('Writing common articles between ...')
    for row in reader:
        if int(row[2]) < 1:
            continue
        print('\t{} and {}'.format(row[0], row[1]))
        idList = p.findall(row[3])
        handle = Entrez.efetch(db="pubmed", id=idList,
                               rettype="medline", retmode="text")
        records = Medline.parse(handle)
        records = list(records)
        for record in records:
            record_file = open(os.path.join(
                output_dir, record.get("PMID", "?") + '.txt'), 'a')
            writer = csv.writer(record_file, delimiter='\t')
            writer.writerow(('PMID: ', record.get("PMID", "?")))
            writer.writerow(('Title: ', record.get("TI", "?")))
            writer.writerow(('Author: ', ', '.join(record.get("AU", "?"))))
            writer.writerow(('Source: ', record.get("SO", "?")))
            writer.writerow(('Abstract: ', record.get("AB", "?")))
            record_file.close()
    assoc_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=False,
                        help="File containing list of genes")
    args = parser.parse_args()
    gene_file = '../analysis/merged.txt'
    gene_synonyms_file = '../data/homo_sapiens_gene_synonyms.txt'
    phenogram_file = '../analysis/merged.txt'
    Entrez.email = "t.fadason@auckland.ac.nz"
    Entrez.api_key = "d3fc89daba8677d8482cc1593c47316e0f08"
    outfile = '../analysis/esearch_metabolites.txt'
    output_dir = '../analysis/downloaded_articles'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    mergedDf, esearch_results = get_gene_synonyms(
        gene_file, gene_synonyms_file)
    phenogramDf, esearch_results = parse_phenogram(
        phenogram_file, esearch_results)
    esearch_results = search_gene_metabolite(
        mergedDf, phenogramDf, esearch_results, outfile)
    write_search_results(esearch_results, outfile)
    download_articles(outfile)
