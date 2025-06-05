import os
import time
import numpy as np
import pandas as pd
from Bio import SeqIO, Entrez
from get_full_sequences_uniprot import get_seq_from_uniprot, get_seq_for_obsolete_id, get_info_from_uniprot, get_url
from utils import sort_id


def filter_hpidb_pairs(output_filepath, filter_by_interaction=False, only_direct_interactions=False, filter_by_confidence=False):
    df = pd.read_table('../data/hpidb2/hpidb2.mitab.txt', sep='\t', header=0,
                       encoding='latin-1').rename(columns={'# protein_xref_1': 'protein_xref_1'})

    df['pair'] = df['protein_xref_1'] + '_' + df['protein_xref_2']

    def get_miscore(x):
        miscore = np.nan
        if not pd.isnull(x):
            scores = x.split('|')
            for score in scores:
                if 'miscore' in score:
                    miscore = float(score.split(':')[-1])
                    break
        return miscore

    df['miscore'] = df['confidence'].apply(get_miscore)

    print('df shape before filtering: ' + str(df.shape))

    # Exclude pairs with detection_method = 'psi-mi:MI:0114(x-ray crystallography)', since they're probably in the PDB
    df = df[~(df['detection_method'] == 'psi-mi:MI:0114(x-ray crystallography)')]
    print('df shape after removing x-ray crystallography pairs: ' + str(df.shape))

    # Exclude pairs from Krogan lab publications (and the same pairs from other papers)
    krogan_lab_papers = ['pubmed:26118995', 'pubmed:30118682', 'pubmed:30550790', 'pubmed:30550789', 'pubmed:22190034',
                         'pubmed:30209081', 'pubmed:25544563', 'pubmed:30833725', 'pubmed:32353859', 'pubmed:33060197',
                         'pubmed:25616068', 'pubmed:37758692', 'pubmed:31527793']
    krogan_lab_papers_df = df[df['pmid'].str.contains('|'.join(krogan_lab_papers), case=False, regex=True)]
    df = df[~df['pair'].isin(krogan_lab_papers_df['pair'].unique().tolist())]
    print('df shape after removing pairs from Krogan lab papers and the same pairs from other papers: ' + str(df.shape))

    if filter_by_interaction:  # only keep direct or physical interactions
        to_exclude = ['psi-mi:MI:0914(association)', 'dip:dip:0914(association)', 'psi-mi:MI:0403(colocalization)',
                   'psi-mi:MI:0796(suppressive genetic interaction defined by inequality)',
                   'psi-mi:MI:0799(additive genetic interaction defined by inequality)',
                   'psi-mi:MI:0794(synthetic genetic interaction defined by inequality)']
        if only_direct_interactions:
            to_exclude.append('psi-mi:MI:0915(physical association)')
        df = df[~(df['interaction_type'].isin(to_exclude))]
    print('df shape after filtering by interaction type: ' + str(df.shape))

    if filter_by_confidence:
        df = df[df['miscore'] >= 0.45] # cutoff recommended by IntAct
    print('df shape after filtering by MIscore: ' + str(df.shape)) # 3514

    # remove duplicates
    print('df shape before removing duplicates: ' + str(df.shape))
    df.sort_values(by=['miscore'], ascending=False, inplace=True)
    df.drop_duplicates(subset=['protein_xref_1', 'protein_xref_2'], keep='first', inplace=True, ignore_index=True)
    print('df shape after removing duplicates: ' + str(df.shape))

    df2 = df[
        ((df['protein_taxid_1'].str.contains('taxid:9606')) & (~df['protein_taxid_2'].str.contains('taxid:9606'))) | (
                (~df['protein_taxid_1'].str.contains('taxid:9606')) & (
            df['protein_taxid_2'].str.contains('taxid:9606')))]

    df2['pathogen_taxid'] = df2['protein_taxid_2'].str.split('(').str[0].str.replace('taxid:', '')
    print('df2 shape: ' + str(df2.shape[0]))

    def get_ids(x):
        xrefs = x.split('|')
        for xref in xrefs:
            if 'uniprotkb:' in xref or 'refseq' in xref:
                protein_id = xref.split(':')[-1]
            elif 'intact:' in xref:
                protein_id = xref.upper()
            else:
                protein_id = ''
        return protein_id

    df2['human_protein_id'] = df2['protein_xref_1'].apply(get_ids)
    df2['pathogen_protein_id'] = df2['protein_xref_2'].apply(get_ids)

    df2.to_csv(output_filepath, index=False)


def get_sequences(hpidb_filepath, output_filepath):
    df = pd.read_csv(hpidb_filepath)
    df_no_nan = df[(~df['human_protein_id'].isna()) & (~df['pathogen_protein_id'].isna())]
    unique_ids = list(set(df_no_nan['human_protein_id'].to_list() + df_no_nan['pathogen_protein_id'].to_list()))
    print(len(unique_ids))
    sequence_dict = {'protein_id': [], 'sequence': []}
    for prot_id in unique_ids: # unique_ids
        print(prot_id)
        sequence_dict['protein_id'].append(prot_id)

        if prot_id.startswith('NP_'):
            seq = get_ncbi_sequence(prot_id, entrez_email_file='../data/hpidb2/entrez_email.txt')
            sequence_dict['sequence'].append(seq)
            time.sleep(1)
        elif prot_id.startswith('INTACT:'):
            seq = get_intact_seqs(prot_id)
            sequence_dict['sequence'] .append(seq)
        elif '-PRO' in prot_id:
            accession = prot_id.split('-')[0]
            chain_id = prot_id.split('-')[-1]
            try:
                uniprot_info = get_info_from_uniprot(accession)
                seq = uniprot_info['results'][0]['sequence']['value']
                for feature in uniprot_info['results'][0]['features']:
                    if feature['type'] == 'Chain':
                        if feature['featureId'] == chain_id:
                            start = int(feature['location']['start']['value'])
                            end = int(feature['location']['end']['value'])
                            break
                subseq = seq[start-1:end]
            except Exception as e:
                print(e)
                subseq = ''
            sequence_dict['sequence'].append(subseq)
            time.sleep(1)
        else:
            try:
                result = get_seq_from_uniprot(prot_id)
                seq = str(result.seq)
            except Exception as e:
                print(e)
                seq = ''
            sequence_dict['sequence'].append(seq)
            time.sleep(1)

    sequence_df = pd.DataFrame(sequence_dict)
    sequence_df.to_csv(output_filepath, index=False)


def get_ncbi_sequence(refseq_id, entrez_email_file):
    # get sequences based on refseq IDs
    with open(entrez_email_file, 'r') as f:
        Entrez.email = f.readlines()[0].strip()

    try:
        handle = Entrez.efetch(db="protein", id=refseq_id, rettype='fasta')
        results = [x.strip() for x in handle.readlines()]
        seq = ''.join(results[1:])
        handle.close()
    except Exception as e:
        print(e)
        seq = ''
    return seq


def get_intact_seqs(intact_id):
    seq_records = SeqIO.to_dict(SeqIO.parse('../data/hpidb2/intact.fasta', 'fasta'))
    try:
        seq = str(seq_records[intact_id.upper()].seq)
    except Exception as e:
        print(e)
        seq = ''
    return seq


def create_input_table(hpidb_filepath, sequences_filepath, output_dir):
    hpidb_df = pd.read_csv(hpidb_filepath)
    hpidb_df_nonan = hpidb_df[(~hpidb_df['human_protein_id'].isna()) & (~hpidb_df['pathogen_protein_id'].isna())]
    print(hpidb_df.shape)
    print(hpidb_df_nonan.shape)
    seq_df = pd.read_csv(sequences_filepath)
    seq_df.set_index(['protein_id'], drop=True, inplace=True)
    seq_dict = seq_df.to_dict()['sequence']

    def get_pair_seq(row):
        human_seq = seq_dict[row['human_protein_id']]
        pathogen_seq = seq_dict[row['pathogen_protein_id']]
        if pd.isnull(human_seq) or pd.isnull(pathogen_seq):
            print('here')
            pair_seq = np.nan
        else:
            pair_seq = f'{human_seq}:{pathogen_seq}'
        return pair_seq

    hpidb_df_nonan['id'] = hpidb_df_nonan['human_protein_id'] + '_' + hpidb_df_nonan['pathogen_protein_id']
    hpidb_df_nonan['sequence'] = hpidb_df_nonan.apply(get_pair_seq, axis=1)
    print(hpidb_df_nonan)
    hpidb_df_nonan.to_csv(os.path.join(output_dir, 'hpidb_selected_pairs_with_seqs.csv'), index=False)
    hpidb_df_nonan.dropna(subset=['sequence'], inplace=True)
    hpidb_yeast = hpidb_df_nonan[hpidb_df_nonan['protein_taxid_2'].str.contains('taxid:559292')]
    hpidb_no_yeast = hpidb_df_nonan[~hpidb_df_nonan['protein_taxid_2'].str.contains('taxid:559292')]
    no_yeast_colabfold_df = hpidb_no_yeast[['id', 'sequence']]
    no_yeast_colabfold_df.to_csv(os.path.join(output_dir, 'other_pathogens', 'colabfold_input.csv'), index=False)
    yeast_colabfold_df = hpidb_yeast[['id', 'sequence']]
    yeast_colabfold_df.to_csv(os.path.join(output_dir, 'scerevisiae', 'colabfold_input.csv'), index=False)


def create_input_table2(hpidb_filepath, sequences_filepath, output_dir):
    hpidb_df = pd.read_csv(hpidb_filepath)
    hpidb_df_nonan = hpidb_df[(~hpidb_df['human_protein_id'].isna()) & (~hpidb_df['pathogen_protein_id'].isna())]
    print(hpidb_df.shape) # 3225
    print(hpidb_df_nonan.shape) # 3246
    seq_df = pd.read_csv(sequences_filepath)
    seq_df.set_index(['protein_id'], drop=True, inplace=True)
    seq_dict = seq_df.to_dict()['sequence']

    def get_pair_seq(row):
        human_seq = seq_dict[row['human_protein_id']]
        pathogen_seq = seq_dict[row['pathogen_protein_id']]
        if pd.isnull(human_seq) or pd.isnull(pathogen_seq):
            pair_seq = np.nan
        else:
            pair_seq = f'{human_seq}:{pathogen_seq}'
        return pair_seq

    hpidb_df_nonan['id'] = hpidb_df_nonan['human_protein_id'] + '_' + hpidb_df_nonan['pathogen_protein_id']
    hpidb_df_nonan['sequence'] = hpidb_df_nonan.apply(get_pair_seq, axis=1)
    hpidb_df_nonan.dropna(subset=['sequence'], inplace=True)
    print(hpidb_df_nonan.shape) # 3225

    # Remove organisms that don't infect humans and are not models commonly used to study related human pathogens
    taxids_to_exclude = [559292, 44689, 644223, 4932, 358769, 11098, 11097, 73475, 284812,
                         12118, 295358, 262719, 11862, 10500, 353796, 10497, 31530, 10389]
    hpidb_selected = hpidb_df_nonan[~hpidb_df_nonan['pathogen_taxid'].isin(taxids_to_exclude)]
    print(hpidb_selected.shape) # 1570 (before filtering like this it was 1632)

    # Remove pairs that were already in the benchmark datasets:
    benchmark_df = pd.read_csv('../data/fullseq_benchmark_dataset_deduplicated.csv')
    benchmark_sorted_ids = benchmark_df['sorted_id'].dropna().unique().tolist()
    hpidb_selected['sorted_id'] = hpidb_selected['id'].apply(sort_id)
    hpidb_selected_no_benchmark = hpidb_selected[~hpidb_selected['sorted_id'].isin(benchmark_sorted_ids)]
    hpidb_selected_no_benchmark.to_csv(os.path.join(output_dir, 'hpidb_selected_pairs_seqs_filtered.csv'), index=False)
    print(hpidb_selected_no_benchmark.shape) # 1536
    colabfold_df = hpidb_selected_no_benchmark[['id', 'sequence']]
    colabfold_df.to_csv(os.path.join(output_dir, 'colabfold_input.csv'), index=False)


def create_input_table3(hpidb_filepath, sequences_filepath, output_dir):
    hpidb_df = pd.read_csv(hpidb_filepath)
    hpidb_df_nonan = hpidb_df[(~hpidb_df['human_protein_id'].isna()) & (~hpidb_df['pathogen_protein_id'].isna())]
    print(hpidb_df.shape) # 3225
    print(hpidb_df_nonan.shape) # 3246
    seq_df = pd.read_csv(sequences_filepath)
    seq_df.set_index(['protein_id'], drop=True, inplace=True)
    seq_dict = seq_df.to_dict()['sequence']

    def get_pair_seq(row):
        human_seq = seq_dict[row['human_protein_id']]
        pathogen_seq = seq_dict[row['pathogen_protein_id']]
        if pd.isnull(human_seq) or pd.isnull(pathogen_seq):
            pair_seq = np.nan
        else:
            pair_seq = f'{human_seq}:{pathogen_seq}'
        return pair_seq

    hpidb_df_nonan['id'] = hpidb_df_nonan['human_protein_id'] + '_' + hpidb_df_nonan['pathogen_protein_id']
    hpidb_df_nonan['sequence'] = hpidb_df_nonan.apply(get_pair_seq, axis=1)
    hpidb_df_nonan.dropna(subset=['sequence'], inplace=True)
    print(hpidb_df_nonan.shape) # 3225

    # Remove organisms that don't infect humans and are not models commonly used to study related human pathogens
    taxids_to_exclude = [559292, 44689, 644223, 4932, 358769, 11098, 11097, 73475, 284812,
                         12118, 295358, 262719, 11862, 10500, 353796, 10497, 31530, 10389]
    hpidb_selected = hpidb_df_nonan[~hpidb_df_nonan['pathogen_taxid'].isin(taxids_to_exclude)]
    print(hpidb_selected.shape) # 1570 (before filtering like this it was 1632)
    hpidb_selected.to_csv(os.path.join(output_dir, 'hpidb_selected_pairs_seqs_filtered.csv'), index=False)
    colabfold_df = hpidb_selected[['id', 'sequence']]
    colabfold_df.to_csv(os.path.join(output_dir, 'colabfold_input.csv'), index=False)


if __name__ == '__main__':
    # filter_hpidb_pairs(output_filepath='../data/hpidb2/hpidb_selected_pairs.csv',
    #                    filter_by_interaction=True, only_direct_interactions=False, filter_by_confidence=True)
    # get_sequences('../data/hpidb2/hpidb_selected_pairs.csv',
    #               '../data/hpidb2/hpidb_selected_pairs_seqs.csv')
    # create_input_table('../data/hpidb2/hpidb_selected_pairs.csv',
    #                    '../data/hpidb2/hpidb_selected_pairs_seqs.csv',
    #                   '../data/hpidb2')
    # create_input_table2('../data/hpidb2/hpidb_selected_pairs.csv',
    #                     '../data/hpidb2/hpidb_selected_pairs_seqs.csv',
    #                     '../data/hpidb2/hpidb_filtered')
    create_input_table3('../data/hpidb2/hpidb_selected_pairs.csv',
                        '../data/hpidb2/hpidb_selected_pairs_seqs.csv',
                        '../data/hpidb2/hpidb_filtered_20250324')
