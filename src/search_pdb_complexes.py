import os
import time
import requests
import json
from bs4 import BeautifulSoup as BS
from urllib.request import urlopen
import re
import numpy as np
import pandas as pd
from Bio import Entrez


def flatten_json(y, exclude=[]):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                if a in exclude:
                    if a == 'rcsb_entity_source_organism':
                        if x[a] is None:
                            out[name + a + '.ncbi_scientific_name'] = np.nan
                            out[name + a + '.ncbi_taxonomy_id'] = np.nan
                        else:
                            out[name + a + '.ncbi_scientific_name'] = ';'.join(
                                [organism['ncbi_scientific_name'] for organism in x[a]])
                            out[name + a + '.ncbi_taxonomy_id'] = ';'.join(
                                [str(organism['ncbi_taxonomy_id']) for organism in x[a]])
                    elif a == 'uniprots':
                        if x[a] is None:
                            out[name + a + '.uniprot_id'] = np.nan
                            out[name + a + '.uniprot_source_organism_taxonomy_id'] = np.nan
                            out[name + a + '.uniprot_source_organism_scientific_name'] = np.nan
                        else:
                            out[name + a + '.uniprot_id'] = ';'.join(
                                [uniprot['rcsb_uniprot_container_identifiers']['uniprot_id'] for uniprot in x[a]])
                            out[name + a + '.uniprot_source_organism_taxonomy_id'] = ';'.join(
                                [str(uniprot['rcsb_uniprot_protein']['source_organism'][
                                         'taxonomy_id']) if 'rcsb_uniprot_protein' in uniprot else '' for uniprot in
                                 x[a]])
                            out[name + a + '.uniprot_source_organism_scientific_name'] = ';'.join(
                                [str(uniprot['rcsb_uniprot_protein']['source_organism'][
                                         'scientific_name']) if 'rcsb_uniprot_protein' in uniprot else '' for
                                 uniprot in x[a]])
                    elif a == 'rcsb_polymer_name_combined':
                        if x[a] is None:
                            out[name + a] = np.nan
                        else:
                            out[name + a] = ';'.join(
                                [polymer_name for polymer_name in x[a]['names']])
                    elif a == 'rcsb_cluster_membership':
                        if x[a] is None:
                            pass
                        else:
                            for seq_cluster in x[a]:
                                out[name + a + '.cluster_id_{}identity'.format(seq_cluster['identity'])] = str(
                                    seq_cluster[
                                        'cluster_id'])
                    else:
                        flatten("; ".join(x[a]), name + a + '.')
                else:
                    flatten(x[a], name + a + '.')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i + 1) + '.')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


def create_search_api_query_assemblies(organism1, organism2, multimer_type, resolution_upper_limit):
    multimer_names_to_num_chains = {'dimer': 2,
                                    'trimer': 3,
                                    'tetramer': 4,
                                    'pentamer': 5,
                                    'hexamer': 6,
                                    'heptamer': 7,
                                    'octamer': 8}

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": organism1,
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": organism2,
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "equals",
                        "value": 2,
                        "attribute": "rcsb_entry_info.polymer_entity_taxonomy_count"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "Hetero {}-mer".format(multimer_names_to_num_chains[multimer_type]),
                        "attribute": "rcsb_struct_symmetry.oligomeric_state"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "equals",
                        "value": multimer_names_to_num_chains[multimer_type],
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_protein"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "equals",
                        "value": 0,
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_nucleic_acid"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "equals",
                        "value": 0,
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_nucleic_acid_hybrid"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "less_or_equal",
                        "value": resolution_upper_limit,
                        "attribute": "rcsb_entry_info.resolution_combined"
                    }
                }
            ]
        },
        "request_options": {
            "results_content_type": [
                "experimental"
            ],
            "return_all_hits": True
        },
        "return_type": "assembly"
    }

    return query


def create_graphql_query_assemblies(assembly_ids):
    format_ids = '[' + ', '.join(['"' + i + '"' for i in assembly_ids]) + ']'
    graphql_query = """{
    assemblies(assembly_ids:""" + format_ids + """){
    rcsb_id
    entry{
      struct {
        title
      }
      rcsb_accession_info{
        initial_release_date
      }
      rcsb_entry_info {
        resolution_combined
        experimental_method
      }
    }
    rcsb_struct_symmetry{
    oligomeric_state
    }
    pdbx_struct_assembly
    {
    oligomeric_details
    rcsb_details
    }
    rcsb_assembly_info{
    polymer_entity_count_protein
    polymer_entity_instance_count_protein
    total_number_interface_residues
    num_protein_interface_entities
    num_heteromeric_interface_entities
    modeled_polymer_monomer_count
    }
    polymer_entity_instances{
      polymer_entity{
        rcsb_polymer_entity{
          pdbx_description
          pdbx_mutation
          rcsb_polymer_name_combined{
            names
          }
        }
        uniprots{
          rcsb_uniprot_container_identifiers{
            uniprot_id
          }
          rcsb_uniprot_protein{
            source_organism{
              taxonomy_id
              scientific_name
            }
          }
        }
        rcsb_entity_source_organism{
          ncbi_taxonomy_id
          ncbi_scientific_name
        }
        entity_poly{
          rcsb_mutation_count
          rcsb_insertion_count
          rcsb_deletion_count
          rcsb_conflict_count
          rcsb_artifact_monomer_count
          rcsb_sample_sequence_length
        }
        rcsb_cluster_membership{
          cluster_id
          identity
        }
      }
    }
    }
    }"""

    return graphql_query


def create_graphql_query_assemblies(assembly_ids):
    format_ids = '[' + ', '.join(['"' + i + '"' for i in assembly_ids]) + ']'
    graphql_query = """{
    assemblies(assembly_ids:""" + format_ids + """){
    rcsb_id
    entry{
      struct {
        title
      }
      rcsb_accession_info{
        initial_release_date
      }
      rcsb_entry_info {
        resolution_combined
        experimental_method
      }
    }
    rcsb_struct_symmetry{
    oligomeric_state
    }
    pdbx_struct_assembly
    {
    oligomeric_details
    rcsb_details
    }
    rcsb_assembly_info{
    polymer_entity_count_protein
    polymer_entity_instance_count_protein
    total_number_interface_residues
    num_protein_interface_entities
    num_heteromeric_interface_entities
    modeled_polymer_monomer_count
    }
    polymer_entity_instances{
      polymer_entity{
        rcsb_polymer_entity{
          pdbx_description
          pdbx_mutation
          rcsb_polymer_name_combined{
            names
          }
        }
        uniprots{
          rcsb_uniprot_container_identifiers{
            uniprot_id
          }
          rcsb_uniprot_protein{
            source_organism{
              taxonomy_id
              scientific_name
            }
          }
        }
        rcsb_entity_source_organism{
          ncbi_taxonomy_id
          ncbi_scientific_name
        }
        entity_poly{
          rcsb_mutation_count
          rcsb_insertion_count
          rcsb_deletion_count
          rcsb_conflict_count
          rcsb_artifact_monomer_count
          rcsb_sample_sequence_length
        }
        rcsb_cluster_membership{
          cluster_id
          identity
        }
      }
    }
    }
    }"""

    return graphql_query


def create_graphql_query_entries(pdb_ids):
    format_pdb_ids = '[' + ', '.join(['"' + i + '"' for i in pdb_ids]) + ']'
    graphql_query = """{
    entries(entry_ids:""" + format_pdb_ids + """) {
    rcsb_id
    struct {
      title
    }
    rcsb_accession_info{
    initial_release_date
    }
    rcsb_entry_info {
      resolution_combined
      experimental_method
    }
    }
  }"""
    return graphql_query


def make_request_search_api(query_dict):
    query_json = json.dumps(query_dict)
    response = requests.get('https://search.rcsb.org/rcsbsearch/v2/query?json={query}'.format(query=query_json))
    if response.status_code == 200:
        return response.json()
    else:
        print('[No data retrieved - %s] %s' % (response.status_code, response.text))
        return {}


def make_request_graphql_api(query):
    graphql_url = 'https://data.rcsb.org/graphql'
    response = requests.post(graphql_url, json={'query': query})
    return response.json()


def get_host_organism_ncbi_taxonomy(tax_id):
    address = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}'.format(tax_id)
    page = urlopen(address).read()
    soup = BS(page, features='html.parser')
    results = soup.body.form.find_all(string=re.compile('Host:.*'), recursive=True)

    if len(results) != 0:
        host = results[0].split(': ')[-1]
    else:
        # print(tax_id)
        host = np.nan

    return host


def remove_nonhuman_pathogens(ids_df):
    cols = [col for col in ids_df.columns.to_list() if 'ncbi_taxonomy_id' in col]
    tax_ids = []
    for col in cols:
        tax_ids.extend(ids_df[col].unique().tolist())
    tax_ids = list(set([id for id in tax_ids if set(id.split(';')) != {'9606'}]))
    # 9606 = Homo sapiens. List comprehension because some entries have '9606; 9606'.
    # Set because some entries can have multiple organisms for a given protein entity (chimeric proteins)

    to_exclude = []
    for tax_id in tax_ids:
        if ';' in tax_id:  # cases where a protein entity has several source organisms
            tax_id_list = tax_id.split(';')
            if len(set(tax_id_list)) != 1:
                # more than 1 source organism for this protein entity, so protein may be chimeric. exclude these cases
                to_exclude.append(tax_id)
                continue
            else:
                tax_id = tax_id_list[0]

        if tax_id in ['223997', '11723',
                      '12643']:  # Murine norovirus 1 (223997), Simian immunodeficiency virus (11723), Ectromelia virus (12643)
            to_exclude.append(tax_id)
            continue

        time.sleep(1)
        host = get_host_organism_ncbi_taxonomy(tax_id)  # get host from NCBI Taxonomy Browser
        if 'human' not in host:
            to_exclude.append(tax_id)
    # print(to_exclude)

    # exclude from ids_df
    filtered_df = ids_df[~ids_df[cols].isin(to_exclude).any(axis=1)]

    return filtered_df


def remove_organisms(ids_df, organisms_to_exclude):
    """Remove rows from ids_df if an organism from organisms_to_exclude is found in that row"""
    filtered_df = ids_df[~ids_df.isin(organisms_to_exclude).any(axis=1)]
    return filtered_df


def reduce_redundancy_uniprot(ids_df):
    uniprot_cols = [col for col in ids_df.columns.to_list() if 'uniprot_id' in col]
    protein_id_cols = []
    for uniprot_col in uniprot_cols:
        names_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.rcsb_polymer_name_combined'.format(
            uniprot_col.split('.')[1])
        description_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.pdbx_description'.format(
            uniprot_col.split('.')[1])
        tax_id_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_entity_source_organism.ncbi_taxonomy_id'.format(
            uniprot_col.split('.')[1])
        protein_id_col = 'polymer_entity_instances.{}.protein_id'.format(uniprot_col.split('.')[1])
        protein_id_cols.append(protein_id_col)
        ids_df[protein_id_col] = ids_df[uniprot_col]
        combined = ids_df[tax_id_col] + '|' + ids_df[names_col]
        ids_df[names_col] = ids_df[names_col].replace(r'^\s*$', np.nan, regex=True).fillna(value=ids_df[
            description_col])  # because some chains don't have an associated rcsb_polymer_name_combined value
        ids_df[protein_id_col] = ids_df[protein_id_col].replace(r'^\s*$', np.nan, regex=True).fillna(
            value=combined)  # because some chains don't have an associated uniprot_id.
        # adding ncbi_taxonomy_id for cases where we don't have the uniprot_id, because different species may have
        # proteins with the same name

    ids_df['complex_id'] = ids_df[protein_id_cols].apply(lambda row: '_'.join(row.values), axis=1)

    complex_ids_to_groups_mapper = {}
    for id in ids_df['complex_id'].unique().tolist():
        inv_id = '_'.join(id.split('_')[::-1])
        if inv_id in complex_ids_to_groups_mapper:
            complex_ids_to_groups_mapper[id] = inv_id
        else:
            complex_ids_to_groups_mapper[id] = id

    ids_df['group'] = ids_df['complex_id'].apply(
        lambda x: complex_ids_to_groups_mapper[x])

    assembly_ids_to_keep = []
    mut_del_ins_cols = [col for col in ids_df.columns.to_list() if
                        ('rcsb_deletion_count' in col) or ('rcsb_insertion_count' in col) or (
                                'rcsb_mutation_count' in col)]
    for group in ids_df['group'].unique().tolist():
        ids_df_subset = ids_df[ids_df['group'] == group]
        if ids_df_subset.shape[0] == 1:
            assembly_ids_to_keep.append(ids_df_subset['assembly_id'].to_list()[0])
        else:
            # exclude mutant entries, unless this results in zero entries:
            ids_df_subset2 = ids_df_subset.loc[(ids_df_subset[mut_del_ins_cols] == 0).all(axis=1)]
            if ids_df_subset2.shape[0] != 0:
                ids_df_subset = ids_df_subset2

            # choose the largest complex using 'rcsb_assembly_info.modeled_polymer_monomer_count'
            largest_structure_df = ids_df_subset[
                ids_df_subset['rcsb_assembly_info.modeled_polymer_monomer_count'] == ids_df_subset[
                    'rcsb_assembly_info.modeled_polymer_monomer_count'].max()]
            if largest_structure_df.shape[0] == 1:
                assembly_ids_to_keep.append(largest_structure_df['assembly_id'].to_list()[0])
            else:  # if more than one entry has the same size, choose the one with the lowest resolution:
                best_resolution = largest_structure_df.loc[
                                  largest_structure_df['entry.rcsb_entry_info.resolution_combined.1'].idxmin(), :]
                assembly_ids_to_keep.append(best_resolution['assembly_id'])

    nonredundant_df = ids_df[ids_df['assembly_id'].isin(assembly_ids_to_keep)]
    # print(len(assembly_ids_to_keep))
    # print(nonredundant_df.shape)
    # print(len(nonredundant_df['complex_id'].to_list()))
    # print(len(nonredundant_df['group'].to_list()))

    return nonredundant_df


def reduce_redundancy_pdbclusters(ids_df, identity_threshold):
    cluster_id_cols = [col for col in sorted(ids_df.columns) if
                       'cluster_id_{}identity'.format(str(identity_threshold)) in col]
    for cluster_id_col in cluster_id_cols:
        uniprot_col = 'polymer_entity_instances.{}.polymer_entity.uniprots.uniprot_id'.format(
            cluster_id_col.split('.')[1])
        polymer_name_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.rcsb_polymer_name_combined'.format(
            cluster_id_col.split('.')[1])
        description_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.pdbx_description'.format(
            uniprot_col.split('.')[1])
        ids_df[cluster_id_col] = ids_df[cluster_id_col].fillna(value=ids_df[uniprot_col]).fillna(
            value=ids_df[polymer_name_col]).fillna(value=ids_df[description_col])
        # small peptides (less than 10 aa) aren't clustered by PDB, so they won't have a cluster_id attributed to them.
        # using uniprot_id instead for these cases, or rcsb_polymer_name_combined for proteins without a Uniprot ID

    ids_df['complex_cluster_ids'] = ids_df[cluster_id_cols].apply(lambda row: '_'.join(row.values), axis=1)
    ids_df['group'] = ids_df[cluster_id_cols].apply(lambda row: '_'.join(sorted(row.values.tolist())), axis=1)
    # ids_df.to_csv('search_results_redundant_with_clusters.csv', index=False)

    assembly_ids_to_keep = []

    for group in ids_df['group'].unique().tolist():
        ids_df_subset = ids_df[ids_df['group'] == group]
        if ids_df_subset.shape[0] == 1:
            assembly_ids_to_keep.append(ids_df_subset['assembly_id'].to_list()[0])
        else:
            # ids_df_subset.to_csv(group + '.csv', index=False)
            # choose the largest complex using 'rcsb_assembly_info.modeled_polymer_monomer_count'
            largest_structure_df = ids_df_subset[
                ids_df_subset['rcsb_assembly_info.modeled_polymer_monomer_count'] == ids_df_subset[
                    'rcsb_assembly_info.modeled_polymer_monomer_count'].max()]
            if largest_structure_df.shape[0] == 1:
                assembly_ids_to_keep.append(largest_structure_df['assembly_id'].to_list()[0])
            else:  # if more than one entry has the same size, choose the one with the lowest resolution:
                best_resolution = largest_structure_df.loc[
                                  largest_structure_df['entry.rcsb_entry_info.resolution_combined.1'].idxmin(), :]
                assembly_ids_to_keep.append(best_resolution['assembly_id'])

    nonredundant_df = ids_df[ids_df['assembly_id'].isin(assembly_ids_to_keep)]

    return nonredundant_df


def reduce_redundancy_pdbclusters2(ids_df, identity_threshold):
    cluster_id_cols = [col for col in sorted(ids_df.columns) if
                       'cluster_id_{}identity'.format(str(identity_threshold)) in col]
    for cluster_id_col in cluster_id_cols:
        uniprot_col = 'polymer_entity_instances.{}.polymer_entity.uniprots.uniprot_id'.format(
            cluster_id_col.split('.')[1])
        ids_df[cluster_id_col] = ids_df[cluster_id_col].fillna(value='NA')
        # small peptides (less than 10 aa) aren't clustered by PDB, so they won't have a cluster_id attributed to them.
        # using uniprot_id instead for these cases, or rcsb_polymer_name_combined for proteins without a Uniprot ID

    ids_df['complex_cluster_ids'] = ids_df[cluster_id_cols].apply(lambda row: '_'.join(row.values), axis=1)
    ids_df['group'] = ids_df[cluster_id_cols].apply(lambda row: '_'.join(sorted(row.values.tolist())), axis=1)

    assembly_ids_to_keep = []
    ids_df_missing_clusters = ids_df[ids_df['group'].str.contains('NA')]
    assembly_ids_to_keep.extend(ids_df_missing_clusters[
                                    'assembly_id'])  # entries that contain small peptides (less than 10 aa) that don't have an associated cluster due to their size
    ids_df_no_missing_clusters = ids_df[~ids_df['group'].str.contains('NA')]

    for group in ids_df_no_missing_clusters['group'].unique().tolist():
        group_df = ids_df_no_missing_clusters[ids_df_no_missing_clusters['group'] == group]
        if group_df.shape[0] == 1:
            assembly_ids_to_keep.append(group_df['assembly_id'].to_list()[0])
        else:
            # group_df.to_csv(group + '.csv', index=False)
            # choose the largest complex using 'rcsb_assembly_info.modeled_polymer_monomer_count'
            largest_structure_df = group_df[
                group_df['rcsb_assembly_info.modeled_polymer_monomer_count'] == group_df[
                    'rcsb_assembly_info.modeled_polymer_monomer_count'].max()]
            if largest_structure_df.shape[0] == 1:
                assembly_ids_to_keep.append(largest_structure_df['assembly_id'].to_list()[0])
            else:  # if more than one entry has the same size, choose the one with the lowest resolution:
                best_resolution = largest_structure_df.loc[
                                  largest_structure_df['entry.rcsb_entry_info.resolution_combined.1'].idxmin(), :]
                assembly_ids_to_keep.append(best_resolution['assembly_id'])

    nonredundant_df = ids_df[ids_df['assembly_id'].isin(assembly_ids_to_keep)]

    return nonredundant_df


def remove_chimeras_and_tags(ids_df):
    cols = [col for col in ids_df.columns if 'uniprot_source_organism_taxonomy_id' in col]

    tax_ids = set()
    for col in cols:
        tax_ids.update(ids_df[col].tolist())

    to_exclude = []
    for tax_id in tax_ids:
        if not pd.isna(tax_id) and ';' in tax_id:
            tax_id_list = tax_id.split(';')
            if len(set(tax_id_list)) != 1:
                # more than 1 source organism for this protein entity, so protein may be chimeric. exclude these cases
                to_exclude.append(tax_id)

    filtered_df = ids_df[~ids_df[cols].isin(to_exclude).any(axis=1)]

    uniprot_id_cols = [col for col in ids_df.columns if 'uniprot_id' in col]
    uniprot_ids = set()
    for uniprot_id_col in uniprot_id_cols:
        uniprot_ids.update(ids_df[uniprot_id_col].tolist())
    to_exclude2 = []
    for uniprot_id in uniprot_ids:
        if not pd.isna(uniprot_id) and ';' in uniprot_id:
            uniprot_id_list = uniprot_id.split(';')
            if len(set(uniprot_id_list)) != 1:
                # more than 1 protein ID for this protein entity, so protein may be chimeric/may have a protein tag added to it. exclude these cases
                to_exclude2.append(uniprot_id)

    filtered_df2 = filtered_df[~filtered_df[uniprot_id_cols].isin(to_exclude2).any(axis=1)]

    filtered_df3 = filtered_df2[~filtered_df2['entry.struct.title'].str.contains('chimeric|chimera', regex=True)]

    return filtered_df3


def remove_missing_organisms(ids_df):
    cols = [col for col in ids_df.columns.to_list() if 'ncbi_taxonomy_id' in col]
    # ids_df[cols] = ids_df[cols].replace(r'^\s*$', np.nan, regex=True)
    filtered_df = ids_df.dropna(axis=0, how='any', subset=cols)
    return filtered_df


def split_data_by_chain_size(ids_df, threshold=30):
    cols = [col for col in ids_df.columns if 'rcsb_sample_sequence_length' in col]
    peptide_protein_df = ids_df[(ids_df[cols] < threshold).any(axis=1)]
    protein_protein_df = ids_df[~(ids_df[cols] < threshold).any(axis=1)]
    # print(peptide_protein_df.shape)
    # print(protein_protein_df.shape)
    # print(set(peptide_protein_df['assembly_id'].to_list()).intersection(set(protein_protein_df['assembly_id'].to_list())))
    return peptide_protein_df, protein_protein_df


def get_taxonomy_lineage(tax_ids):
    Entrez.email = 'dbaptista@igc.gulbenkian.pt'
    Entrez.sleep_between_tries = 30
    Entrez.max_tries = 5
    tax_ids_new_list = []
    multiple_tax_ids_mapper = {}
    for tax_id in tax_ids:
        if isinstance(tax_id, float):
            if np.isnan(tax_id):
                continue
            else:
                tax_id = str(int(tax_id))

        if len(tax_id.split(';')) > 1:
            new_tax_id = tax_id.split(';')[0]
            tax_ids_new_list.append(new_tax_id)
            multiple_tax_ids_mapper[tax_id] = new_tax_id
        else:
            tax_ids_new_list.append(tax_id)

    search = Entrez.efetch(id=tax_ids_new_list, db='taxonomy', retmode='xml')
    records = Entrez.read(search)
    lineages = {}
    for record in records:
        lineages[record['TaxId']] = [x['TaxId'] for x in record['LineageEx']]
    for multiple_tax_id in multiple_tax_ids_mapper:
        lineages[multiple_tax_id] = lineages[multiple_tax_ids_mapper[multiple_tax_id]]

    missing_tax_ids = set(tax_ids_new_list).difference(set(lineages.keys()))  # IDs that redirect to different tax IDs
    for missing_tax_id in missing_tax_ids:
        search = Entrez.efetch(id=missing_tax_id, db='taxonomy', retmode='xml')
        record = Entrez.read(search)[0]
        lineages[missing_tax_id] = [x['TaxId'] for x in record['LineageEx']]

    return lineages


def check_organisms(ids_df):
    ncbi_tax_cols = sorted([col for col in ids_df.columns if 'ncbi_taxonomy_id' in col])
    uniprot_tax_cols = sorted([col for col in ids_df.columns if 'uniprot_source_organism_taxonomy_id' in col])

    all_tax_ids = set()
    for ncbi_tax_col, uniprot_tax_col in zip(ncbi_tax_cols, uniprot_tax_cols):
        all_tax_ids.update(ids_df[ncbi_tax_col])
        all_tax_ids.update(ids_df[uniprot_tax_col])

    tax_lineages = get_taxonomy_lineage(all_tax_ids)
    to_exclude = set()
    for ncbi_tax_col, uniprot_tax_col in zip(ncbi_tax_cols, uniprot_tax_cols):
        lineage_col = 'polymer_entity_instances.{}.polymer_entity.taxonomy_lineage'.format(ncbi_tax_col.split('.')[1])
        uniprot_lineage_col = 'polymer_entity_instances.{}.polymer_entity.uniprots.uniprot_taxonomy_lineage'.format(
            uniprot_tax_col.split('.')[1])
        ids_df[lineage_col] = ids_df[ncbi_tax_col].apply(
            lambda x: tax_lineages[x] if not pd.isna(x) else np.nan)
        ids_df[uniprot_lineage_col] = ids_df[uniprot_tax_col].apply(
            lambda x: tax_lineages[x] if not pd.isna(x) else np.nan)

        for i, row in ids_df.iterrows():
            if (not pd.isna(row[uniprot_tax_col])) and (
                    row[ncbi_tax_col].split(';')[0] != row[uniprot_tax_col].split(';')[0]):
                if (not row[ncbi_tax_col].split(';')[0] in row[uniprot_lineage_col]) and (
                        not row[uniprot_tax_col].split(';')[0] in row[lineage_col]):
                    if not (row[lineage_col] == row[uniprot_lineage_col]):
                        ind = min(len(row[lineage_col]), len(row[uniprot_lineage_col]))
                        if not (row[lineage_col][0:ind] == row[uniprot_lineage_col][0:ind]) and (
                                not row[lineage_col][0:ind - 1] == row[uniprot_lineage_col][0:ind - 1]) and (
                                not row[lineage_col][0:ind - 2] == row[uniprot_lineage_col][0:ind - 2]):
                            to_exclude.add(row['assembly_id'])
    filtered_df = ids_df[~ids_df['assembly_id'].isin(to_exclude)]

    return filtered_df


def get_pdb_ids(host, pathogen, multimer_type='dimer', resolution_upper_limit=3.0, only_keep_first_bioassembly=True,
                exclude_nonhuman_pathogens=False, organisms_to_exclude=None, exclude_software_defined=False,
                min_date=None, output_dir=None):
    query = create_search_api_query_assemblies(organism1=host, organism2=pathogen, multimer_type=multimer_type,
                                               resolution_upper_limit=resolution_upper_limit)
    results = make_request_search_api(query_dict=query)

    df = pd.DataFrame(data=results['result_set'])
    df['pdb_id'] = df['identifier'].apply(lambda x: x.split('-')[0])
    df.rename({'identifier': 'assembly_id'}, axis=1, inplace=True)
    df.drop('score', axis=1, inplace=True)

    if only_keep_first_bioassembly:
        df = df.drop_duplicates(subset='pdb_id', keep='first')

    # get more info on the assemblies using the GraphQL API:
    assemblies_query = create_graphql_query_assemblies(df['assembly_id'].tolist())
    graphql_assemblies_results = make_request_graphql_api(assemblies_query)
    assemblies_info_df = pd.DataFrame(
        [flatten_json(x, exclude=['uniprots', 'rcsb_entity_source_organism', 'rcsb_polymer_name_combined',
                                  'rcsb_cluster_membership']) for x in
         graphql_assemblies_results['data']['assemblies']])

    results_df = df.merge(assemblies_info_df, how='left', left_on='assembly_id', right_on='rcsb_id').drop(['rcsb_id'],
                                                                                                          axis=1)
    assert results_df['assembly_id'].tolist() == df['assembly_id'].tolist()

    if output_dir is not None:
        results_df.to_csv(os.path.join(output_dir, 'full_search_results.csv'), index=False)

    # remove entries with chimeric proteins and experimental tags:
    results_df = remove_chimeras_and_tags(results_df)
    results_df.to_csv(os.path.join(output_dir, 'results_without_chimeras_tags.csv'), index=False)

    # remove entries with polymer entities without an associated source organism (rcsb ncbi_taxonomy_id):
    results_df = remove_missing_organisms(results_df)
    results_df.to_csv(os.path.join(output_dir, 'results_without_prots_with_missing_organisms.csv'), index=False)

    # check organisms and remove entries where the source organisms in RCSB PDB differ from the source organisms of the
    # associated Uniprot IDs
    results_df = check_organisms(results_df)
    results_df.to_csv(os.path.join(output_dir, 'results_without_source_organism_conflicts.csv'), index=False)

    if exclude_nonhuman_pathogens:
        results_df = remove_nonhuman_pathogens(results_df)

    if organisms_to_exclude is not None:
        results_df = remove_organisms(results_df, organisms_to_exclude)

    if exclude_software_defined:
        results_df = results_df[results_df['pdbx_struct_assembly.rcsb_details'] != 'software_defined_assembly']

    if min_date is not None:
        results_df = results_df[
            results_df['entry.rcsb_accession_info.initial_release_date'] >= '{}T:00:00:00Z'.format(min_date)]

    results_df_90identity = reduce_redundancy_pdbclusters2(results_df, identity_threshold=90)
    results_df_70identity = reduce_redundancy_pdbclusters2(results_df, identity_threshold=70)
    assert len(set(results_df_70identity['assembly_id'].to_list()).difference(
        set(results_df_90identity['assembly_id'].tolist()))) == 0

    if output_dir is not None:
        save_results(df=results_df, prefix='filtered', output_dir=output_dir)
        save_results(df=results_df_90identity, prefix='90identity', output_dir=output_dir)
        save_results(df=results_df_70identity, prefix='70identity', output_dir=output_dir)

    return results_df


def save_results(df, output_dir, prefix=''):
    if prefix:
        prefix = '{}_'.format(prefix)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    df.to_csv(os.path.join(output_dir, '{}search_results.csv'.format(prefix)), index=False)

    with open(os.path.join(output_dir, '{}selected_pdb_ids.txt'.format(prefix)), 'w') as f:
        f.write(','.join(df['pdb_id'].unique().tolist()))

    with open(os.path.join(output_dir, '{}selected_assembly_ids.txt'.format(prefix)), 'w') as f:
        f.write(','.join(df['assembly_id'].unique().tolist()))

# if __name__ == '__main__':
#     get_pdb_ids(host='Mammalia', pathogen='Viruses', multimer_type='dimer',
#                 resolution_upper_limit=3.0, only_keep_first_bioassembly=True,
#                 output_dir='../data/virus_mammalia_dimers')
#     # get_pdb_ids(host='Mammalia', pathogen='Bacteria', multimer_type='dimer',
#     #                                resolution_upper_limit=3.0, only_keep_first_bioassembly=True,
#     #                                output_dir='../data/bacteria_mammalia_dimers')
