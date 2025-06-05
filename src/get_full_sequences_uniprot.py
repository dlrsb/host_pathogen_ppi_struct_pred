import os
import copy
import time
import argparse
import pickle
from io import StringIO
import pandas as pd
from Bio import SeqIO
from fuzzywuzzy import fuzz
from search_pdb_complexes import make_request_graphql_api, flatten_json
from utils import load_structure, get_url
from get_uniprot_info import get_info_from_uniprot


def create_graphql_query_chainids(assembly_ids):
    format_ids = '[' + ', '.join(['"' + i + '"' for i in assembly_ids]) + ']'
    graphql_query = """{
    assemblies(assembly_ids:""" + format_ids + """){
    rcsb_id
    polymer_entity_instances{
      polymer_entity{
        rcsb_polymer_entity{
          pdbx_description
          rcsb_polymer_name_combined{
            names
          }
        }
        entity_poly{
          pdbx_strand_id
        }
        uniprots{
          rcsb_uniprot_container_identifiers{
            uniprot_id
          }
        }
      }
    }
    }
    }"""

    return graphql_query


def get_seq_for_obsolete_id(uniprot_id):
    url = f'https://rest.uniprot.org/unisave/{uniprot_id}?format=fasta'
    response = get_url(url)
    seq_record = list(SeqIO.parse(StringIO(response.text), 'fasta'))[0] # the first sequence in the FASTA file is the sequence from the last Uniprot release that included the now obseolete ID
    return seq_record


def get_seq_from_uniprot(uniprot_id):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta'
    response = get_url(url)
    if len(response.text) == 0: # if response.text is empty, download the FASTA from the last version that included the obsolete Uniprot ID
        seq_record = get_seq_for_obsolete_id(uniprot_id)
    else:
        seq_record = list(SeqIO.parse(StringIO(response.text), 'fasta'))[0]

    return seq_record


def check_if_is_polyprotein(uniprot_info):
    names = []
    if 'recommendedName' in uniprot_info['proteinDescription']:
        names.append(uniprot_info['proteinDescription']['recommendedName']['fullName']['value'])
        if 'shortNames' in uniprot_info['proteinDescription']['recommendedName']:
            for name in uniprot_info['proteinDescription']['recommendedName']['shortNames']:
                names.append(name['value'])
    if 'submissionNames' in uniprot_info['proteinDescription']:
        for name in uniprot_info['proteinDescription']['submissionNames']:
            names.append(name['fullName']['value'])
    if 'alternativeNames' in uniprot_info['proteinDescription']:
        for name in uniprot_info['proteinDescription']['alternativeNames']:
            names.append(name['fullName']['value'])

    if any('polyprotein' in x.lower() for x in names):
        return True
    else:
        return False


def get_chain_info(data_df):
    query = create_graphql_query_chainids(data_df['assembly_id'].to_list())
    graphql_results = make_request_graphql_api(query)
    chain_info_df = pd.DataFrame([flatten_json(x, exclude=['uniprots', 'rcsb_polymer_name_combined']) for x in
                                  graphql_results['data']['assemblies']])
    chain_info_df.rename(mapper={'rcsb_id': 'assembly_id'}, axis=1, inplace=True)
    return chain_info_df


def get_chain_id(prot_ids, chain_id_pdb_info, uniprot_cols, native_chains):
    for uniprot_col in uniprot_cols:
        entity_id = uniprot_col.split('.')[1]
        name_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.rcsb_polymer_name_combined'.format(
            entity_id)
        description_col = 'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.pdbx_description'.format(
            entity_id)
        if chain_id_pdb_info[uniprot_col].values in prot_ids or chain_id_pdb_info[name_col].values in prot_ids or \
                chain_id_pdb_info[description_col].values in prot_ids:
            chain_ids_pdb = chain_id_pdb_info[
                'polymer_entity_instances.{}.polymer_entity.entity_poly.pdbx_strand_id'.format(entity_id)]
            chain_ids_pdb_list = [x[0] for x in chain_ids_pdb.to_list()[0].split(
                ',')]  # doing this because some PDB entries have chain IDs like 'EEE', 'FFF', etc.
            # print('Entity ' + str(entity_id))
            # print('Native PDB file: ' + str(native_chains))
            # print('PDB info: ' + str(chain_ids_pdb_list))
            chain_id = set(chain_ids_pdb_list).intersection(set(native_chains)).pop()
            # print('Chain ID: ' + chain_id)
            break

    return chain_id


def get_full_sequences_from_uniprot(dataset_filepath, native_pdb_dir, fasta_output_dir=None):
    df = pd.read_csv(dataset_filepath)
    chain_info_df = get_chain_info(df)
    # df = df.merge(chain_info_df, how='left', on='assembly_id')
    uniprot_cols = [col for col in df.columns.to_list() if 'uniprot_id' in col]

    polyprotein_entries = []

    for i, row in df.iterrows():
        entry_has_polyprotein = False
        id = row['assembly_id']
        # print(id)
        if os.path.exists(os.path.join(native_pdb_dir, '{}_fixedmodel.pdb'.format(id))):
            native_pdb_filepath = os.path.join(native_pdb_dir, '{}_fixedmodel.pdb'.format(id))
        else:
            native_pdb_filepath = os.path.join(native_pdb_dir, '{}.pdb'.format(id))
        # records = list(SeqIO.parse(native_pdb_filepath, 'pdb-atom'))
        # native_chains = [records[i].id.split(':')[-1] for i in
        #           range(len(records))]  # get the chain IDs that are present in this biological assembly file
        struct = load_structure(native_pdb_filepath)
        native_chains = [chain.id for chain in struct]
        seq_records = []
        chain_to_id_mapper = {}
        for col in uniprot_cols:
            is_polyprotein = False
            entity_id = col.split('.')[1]
            chain_id_pdb_info = chain_info_df[chain_info_df['assembly_id'] == id]
            all_protein_ids = [row[x] for x in [col,
                                                'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.rcsb_polymer_name_combined'.format(
                                                    entity_id),
                                                'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.pdbx_description'.format(
                                                    entity_id)] if not pd.isnull(row[x])]
            chain_id = get_chain_id(prot_ids=all_protein_ids, chain_id_pdb_info=chain_id_pdb_info,
                                    uniprot_cols=uniprot_cols, native_chains=native_chains)
            # chain_id_pdb = row['polymer_entity_instances.{}.polymer_entity.entity_poly.pdbx_strand_id'.format(entity_id)]
            # chain_id = set(chain_id_pdb.split(',')).intersection(set(native_chains)).pop()

            # get full sequences if protein has Uniprot ID, else extract sequence from the native PDB file/FASTA file
            # extracted from the PDB file
            if pd.isnull(row[col]):
                # protein_id = row[
                #     'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.rcsb_polymer_name_combined'.format(
                #         entity_id)]
                protein_id = all_protein_ids[0]
                protein_id = protein_id.replace("'", "").replace("/", " ").replace("(", " ").replace(")", " ").replace(",", " ")
                protein_id = ''.join(protein_id.split(' '))
                if pd.isnull(protein_id):
                    protein_id = row[
                        'polymer_entity_instances.{}.polymer_entity.rcsb_polymer_entity.pdbx_description'.format(
                            entity_id)]

                for record in SeqIO.parse(native_pdb_filepath,
                                          "pdb-seqres"):  # get seqs from the SEQRES records in the PDB file
                    if record.annotations['chain'] == chain_id:
                        chain_seq_record = copy.deepcopy(record)
                        chain_seq_record.id = protein_id
                        chain_seq_record.description = ''
                        chain_seq_record_native_chain_id = copy.deepcopy(record)
                        chain_seq_record_native_chain_id.description = '{}_{}'.format(protein_id, id)
                        chain_to_id_mapper[chain_id] = '{}_{}'.format(protein_id, id)
                        single_prot_fasta_filepath = os.path.join(fasta_output_dir, 'individual_proteins',
                                                                  '{}_{}.fasta'.format(protein_id, id))
            else:
                protein_id = row[col]
                chain_seq_record = get_seq_from_uniprot(protein_id)
                chain_seq_record.id = protein_id
                uniprot_info = get_info_from_uniprot(protein_id)
                is_polyprotein = check_if_is_polyprotein(uniprot_info)
                if is_polyprotein:
                    entry_has_polyprotein = True
                    polyprotein_entries.append(id)

                chain_seq_record_native_chain_id = copy.deepcopy(chain_seq_record)
                chain_seq_record_native_chain_id.id = chain_id
                chain_seq_record_native_chain_id.description = protein_id
                chain_to_id_mapper[chain_id] = protein_id
                single_prot_fasta_filepath = os.path.join(fasta_output_dir, 'individual_proteins',
                                                          '{}.fasta'.format(protein_id))
            if not is_polyprotein:
                with open(single_prot_fasta_filepath, 'w') as output_handle:
                    SeqIO.write([chain_seq_record], output_handle, 'fasta-2line')

                seq_records.append(chain_seq_record_native_chain_id)

                time.sleep(1)

        if not entry_has_polyprotein:
            multimer_fasta_path = os.path.join(fasta_output_dir, 'protein_pairs', '{}.fasta'.format(id))
            # sort chains according to order in native file:
            seq_records_sorted = [sr for native_chain_id in native_chains for sr in seq_records if sr.id == native_chain_id]
            with open(multimer_fasta_path, 'w') as output_handle:
                SeqIO.write(seq_records_sorted, output_handle, 'fasta-2line')

            with open(os.path.join(fasta_output_dir, 'chain_maps', '{}_id_to_chain_map.pkl'.format(id)), 'wb') as f:
                pickle.dump(chain_to_id_mapper, f)

    with open(os.path.join(fasta_output_dir, 'entries_with_polyproteins.txt'), 'w') as f:
        f.writelines([line + '\n' for line in polyprotein_entries])


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Get full sequences for Uniprot IDs in a data table and save as FASTA files')
    arg_parser.add_argument(
        '-i',
        dest='dataset_filepath',
        default=None,
        help='Path to the dataset file '
    )
    arg_parser.add_argument(
        '-n',
        dest='native_pdb_dir',
        default=None,
        help='Path to the directory containing the native PDB files'
    )
    arg_parser.add_argument(
        '-o',
        dest='fasta_output_dir',
        default=None,
        help='Path to the output directory where the FASTA files will be saved'
    )
    args = arg_parser.parse_args()
    get_full_sequences_from_uniprot(**vars(args))