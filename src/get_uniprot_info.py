import sys
import time
import requests
import pandas as pd
from utils import get_url


def get_info_from_uniprot(uniprot_id):
    url = 'https://rest.uniprot.org/uniprotkb/accessions?accessions={}'.format(uniprot_id)
    response = get_url(url)
    result = response.json()
    return result


def get_gene_and_protein_names(filepath, column):
    df = pd.read_csv(filepath)
    protein_ids = df[column].unique().tolist()
    results_dict = {'id': [], 'gene': [], 'protein': [], 'organism': []}
    for prot_id in protein_ids:
        print(prot_id)
        try:
            uniprot_result = get_info_from_uniprot(prot_id)
        except requests.exceptions.HTTPError as e:
            print(e)
            print(prot_id)
        else:
            if len(uniprot_result) != 0:
                results_dict['id'].append(prot_id)
                genes = uniprot_result['results'][0]['genes']
                gene_name = ';'.join([genes[i]['geneName']['value'] if 'geneName' in genes[i] else genes[i]['orfNames'][0]['value'] for i in range(len(genes))])
                protein_name = uniprot_result['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                organism = uniprot_result['results'][0]['organism']['scientificName']
                results_dict['gene'].append(gene_name)
                results_dict['protein'].append(protein_name)
                results_dict['organism'].append(organism)
                print(gene_name)
                print(protein_name)
                print(organism)
            time.sleep(1)

    uniprot_df = pd.DataFrame(results_dict)
    df_merged = df.merge(uniprot_df, how='left', left_on=column, right_on='id')
    return df_merged
