import copy
import glob
import os
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def create_individual_fastas_hpidb(hpidb_dataset_path, output_dir):
    hpidb_df = pd.read_csv(hpidb_dataset_path)
    hpidb_df['pathogen_protein_id'] = hpidb_df['pathogen_protein_id'].str.replace('INTACT:', 'INTACT-').str.replace(
        'PRO_', 'PRO-').str.replace('NP_', 'NP-')
    hpidb_df['human_protein_id'] = hpidb_df['human_protein_id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_',
                                                                                                              'PRO-').str.replace(
        'NP_', 'NP-')
    for i, r in hpidb_df.iterrows():
        prot1_id = r['human_protein_id']
        prot2_id = r['pathogen_protein_id']
        prot1_seq, prot2_seq = r['sequence'].split(':')
        record1 = SeqRecord(Seq(prot1_seq), id=prot1_id, name='', description='')
        record2 = SeqRecord(Seq(prot2_seq), id=prot2_id, name='', description='')
        prot1_output_filepath = os.path.join(output_dir, f'{prot1_id}.fasta')
        prot2_output_filepath = os.path.join(output_dir, f'{prot2_id}.fasta')
        if not os.path.exists(prot1_output_filepath):
            SeqIO.write(record1, prot1_output_filepath, 'fasta-2line')
        if not os.path.exists(prot2_output_filepath):
            SeqIO.write(record2, prot2_output_filepath, 'fasta-2line')


def create_fasta_for_mmseqs(virus_benchmark_seqs_path, bacteria_benchmark_seqs_path, apms_seqs_path,
                            sars_mers_sars2_path, hpidb_seqs_path, output_filepath, minimum_length=10):
    all_files = []

    for pth in [virus_benchmark_seqs_path, bacteria_benchmark_seqs_path, apms_seqs_path, sars_mers_sars2_path,
                hpidb_seqs_path]:
        filepaths = [os.path.join(pth, x) for x in os.listdir(pth)]
        all_files.extend(filepaths)

    seen = set()
    records = []
    for fasta in all_files:
        record = SeqIO.read(fasta, 'fasta')
        new_id = os.path.splitext(os.path.split(fasta)[-1])[
            0]  # because for the benchmark proteins sometimes the filename doesn't match the protein ID
        record.id = new_id
        print(record.id)
        record.name = ''
        record.description = ''
        if record.id not in seen and len(record.seq) > minimum_length:
            seen.add(record.id)
            records.append(record)

    SeqIO.write(records, output_filepath, 'fasta-2line')
    print(len(all_files))
    print(len(seen))
    print(len(records))


def get_pair_group(r, cluster_col1, cluster_col2):
    group = '_'.join(sorted([r[cluster_col1], r[cluster_col2]]))
    return group


def add_clusters(pairs_df, col1, col2, clusters_df):
    pairs_clusters_df1 = pairs_df.merge(clusters_df, how='left', left_on=col1, right_on='cluster-member')
    pairs_clusters_df1.rename(columns={'cluster-representative': f'{col1}_cluster'}, inplace=True)
    pairs_clusters_df1.drop(columns=['cluster-member'], inplace=True)
    pairs_clusters_df2 = pairs_clusters_df1.merge(clusters_df, how='left', left_on=col2, right_on='cluster-member')
    pairs_clusters_df2.rename(columns={'cluster-representative': f'{col2}_cluster'}, inplace=True)
    pairs_clusters_df2.drop(columns=['cluster-member'], inplace=True)
    pairs_clusters_df2[[f'{col1}_cluster', f'{col2}_cluster']] = pairs_clusters_df2[
        [f'{col1}_cluster', f'{col2}_cluster']].fillna(value='NA')
    pairs_clusters_df2['pair_cluster'] = pairs_clusters_df2.apply(get_pair_group,
                                                                  args=(f'{col1}_cluster', f'{col2}_cluster'), axis=1)
    return pairs_clusters_df2


def add_benchmark_clusters(benchmark_pair_seqs_dir, seq_clusters_filepath, output_filepath):
    benchmark_pair_files = os.listdir(benchmark_pair_seqs_dir)
    pairs_dict = {'pair_id': [], 'protein1': [], 'protein2': []}
    for pair_file in benchmark_pair_files:
        pairs_dict['pair_id'].append(os.path.splitext(pair_file)[0])
        records = list(SeqIO.parse(os.path.join(benchmark_pair_seqs_dir, pair_file), 'fasta'))
        for i in range(len(records)):
            pairs_dict[f'protein{i + 1}'].append(records[i].description.split(' ')[-1])
    pairs_df = pd.DataFrame(data=pairs_dict)
    # Add clusters
    clusters_df = pd.read_table(seq_clusters_filepath, sep='\t', header=None,
                                names=['cluster-representative', 'cluster-member'])
    pairs_clusters_df = add_clusters(pairs_df=pairs_df, col1='protein1', col2='protein2', clusters_df=clusters_df)
    pairs_clusters_df.to_csv(output_filepath, index=False)


def remove_duplicates_fullseq_benchmark_metadata(bacteria_metadata_filepath, bacteria_clusters_filepath,
                                                 virus_metadata_filepath, virus_clusters_filepath,
                                                 output_filepath):
    bacteria_metadata_df = pd.read_csv(bacteria_metadata_filepath)
    bacteria_clusters_df = pd.read_csv(bacteria_clusters_filepath)
    virus_metadata_df = pd.read_csv(virus_metadata_filepath)
    virus_clusters_df = pd.read_csv(virus_clusters_filepath)
    bacteria_metadata_df['dataset'] = 'bacteria'
    virus_metadata_df['dataset'] = 'virus'
    bacteria_metadata_df = bacteria_metadata_df.merge(bacteria_clusters_df, how='left', left_on='assembly_id',
                                                      right_on='pair_id')
    virus_metadata_df = virus_metadata_df.merge(virus_clusters_df, how='left', left_on='assembly_id',
                                                right_on='pair_id')

    both_df = pd.concat([bacteria_metadata_df, virus_metadata_df], axis=0, ignore_index=True)
    print(both_df.shape)
    print(both_df.columns)
    both_df.sort_values(by='rcsb_assembly_info.modeled_polymer_monomer_count', ascending=False, inplace=True,
                        ignore_index=True)  # to keep the PDB entry with the most modeled residues when deduplicating
    both_df_deduplicated = both_df.drop_duplicates(subset=['pair_cluster'], keep='first', ignore_index=True)
    print(both_df_deduplicated.shape)
    both_df_deduplicated.to_csv(output_filepath, index=False)


def remove_duplicates_benchmark_results(results_filepath, deduplicated_metadata_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    metadata_df = pd.read_csv(deduplicated_metadata_filepath)
    results_df_deduplicated = copy.deepcopy(results_df[results_df['name'].isin(metadata_df['assembly_id'].to_list())])
    results_df_deduplicated.to_csv(output_filepath, index=False)


def remove_duplicates_hpidb_metadata(hpidb_metadata_filepath, seq_id_70_clusters_filepath, benchmark_metadata_filepath,
                                     seq_id_95_clusters_filepath, apms_filepath, output_filepath):
    hpidb_metadata_df = pd.read_csv(hpidb_metadata_filepath)
    hpidb_metadata_df['id'] = hpidb_metadata_df['id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_',
                                                                                                    'PRO-').str.replace(
        'NP_', 'NP-')
    hpidb_metadata_df['pathogen_protein_id'] = hpidb_metadata_df['pathogen_protein_id'].str.replace('INTACT:',
                                                                                                    'INTACT-').str.replace(
        'PRO_', 'PRO-').str.replace('NP_', 'NP-')
    hpidb_metadata_df['human_protein_id'] = hpidb_metadata_df['human_protein_id'].str.replace('INTACT:',
                                                                                              'INTACT-').str.replace(
        'PRO_', 'PRO-').str.replace('NP_', 'NP-')

    print(hpidb_metadata_df.shape)

    hpidb_metadata_df['protein1'] = hpidb_metadata_df['human_protein_id']
    hpidb_metadata_df['protein2'] = hpidb_metadata_df['pathogen_protein_id']

    # Remove pairs that already exist in the benchmark dataset (based on 70% seq id clusters):
    clusters_seqid70_df = pd.read_table(seq_id_70_clusters_filepath, sep='\t', header=None,
                                        names=['cluster-representative', 'cluster-member'])
    benchmark_metadata_df = pd.read_csv(benchmark_metadata_filepath)
    hpidb_metadata_df = add_clusters(pairs_df=hpidb_metadata_df, col1='protein1', col2='protein2',
                                     clusters_df=clusters_seqid70_df)
    hpidb_metadata_df.rename(columns={'protein1_cluster': 'protein1_cluster_seqid70',
                                      'protein2_cluster': 'protein2_cluster_seqid70',
                                      'pair_cluster': 'pair_cluster_seqid70'}, inplace=True)
    hpidb_metadata_df_deduplicated1 = hpidb_metadata_df[
        ~hpidb_metadata_df['pair_cluster_seqid70'].isin(benchmark_metadata_df['pair_cluster'].to_list())]
    print(hpidb_metadata_df_deduplicated1.shape)

    # Remove pairs that already exist in the AP-MS datasets (based on 95% seq id clusters):
    apms_df = pd.read_csv(apms_filepath)
    clusters_seqid95_df = pd.read_table(seq_id_95_clusters_filepath, sep='\t', header=None,
                                        names=['cluster-representative', 'cluster-member'])
    hpidb_metadata_df_deduplicated1 = add_clusters(pairs_df=hpidb_metadata_df_deduplicated1, col1='protein1',
                                                   col2='protein2', clusters_df=clusters_seqid95_df)
    hpidb_metadata_df_deduplicated1.rename(columns={'protein1_cluster': 'protein1_cluster_seqid95',
                                                    'protein2_cluster': 'protein2_cluster_seqid95',
                                                    'pair_cluster': 'pair_cluster_seqid95'}, inplace=True)
    hpidb_metadata_df_deduplicated2 = hpidb_metadata_df_deduplicated1[
        ~hpidb_metadata_df_deduplicated1['pair_cluster_seqid95'].isin(apms_df['pair_cluster_seqid95'].to_list())]
    print(hpidb_metadata_df_deduplicated2.shape)

    # Remove pairs that are duplicates within the HPIDB dataset (based on 95% sequence identity)
    hpidb_metadata_df_deduplicated3 = hpidb_metadata_df_deduplicated2.drop_duplicates(subset=['pair_cluster_seqid95'],
                                                                                      keep='first')
    print(hpidb_metadata_df_deduplicated3.shape)

    hpidb_metadata_df_deduplicated3.to_csv(output_filepath, index=False)


def remove_duplicates_hpidb_results(results_filepath, deduplicated_metadata_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    print(results_df.shape)
    results_df['name'] = results_df['name'].str.replace('INTACT_', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace(
        'NP_', 'NP-')
    metadata_df = pd.read_csv(deduplicated_metadata_filepath)
    results_df_deduplicated = copy.deepcopy(results_df[results_df['name'].isin(metadata_df['id'].to_list())])
    print(results_df_deduplicated.shape)
    results_df_deduplicated.to_csv(output_filepath, index=False)


def apms_input_pairs(output_filepath):
    prot_pairs1 = [x.split('.')[0] for x in
                   os.listdir('../data/krogan_lab_host_pathogen_data/fastafiles/protein_pairs')]
    prot_pairs2 = [x.split('.')[0] for x in
                   os.listdir('../data/krogan_lab_host_pathogen_data/sars_mers_sars2_human_ppi/protein_pairs')]
    all_prot_pairs = prot_pairs1 + prot_pairs2
    df = pd.DataFrame(data={'protein_pair': all_prot_pairs})
    df[['protein1', 'protein2']] = df['protein_pair'].str.split(pat='_', expand=True)

    # Exclude pairs from unpublihsed EV71 dataset
    df_filtered = df[~df['protein_pair'].str.startswith('EV71')]

    # Exclude pairs from unpublished IAV-HTBE dataset that are also not in the other IAV datasets
    iav_htbe_input_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-HTBE.csv')
    iav_a549_input_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-A549.csv')
    iav_thp1_input_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-THP1.csv')
    iav_hek293t_input_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-HEK293T.csv')
    iav_published_df = pd.concat([iav_a549_input_df, iav_thp1_input_df, iav_hek293t_input_df], axis=0,
                                 ignore_index=True)
    iav_pairs_exclude = set(iav_htbe_input_df['protein_pair'].to_list()).difference(
        set(iav_published_df['protein_pair'].to_list()))
    df_filtered2 = df_filtered[~df_filtered['protein_pair'].isin(iav_pairs_exclude)]
    print(df.shape)
    print(df_filtered.shape)
    print(df_filtered2.shape)
    df_filtered2.to_csv(output_filepath, index=False)


def remove_duplicates_apms_input_pairs(input_pairs_filepath, seq_id_70_clusters_filepath, benchmark_metadata_filepath,
                                       seq_id_95_clusters_filepath, output_filepath):
    df = pd.read_csv(input_pairs_filepath)

    # Remove pairs that already exist in the benchmark dataset (based on 70% seq id clusters):
    print(df.shape)
    clusters_seqid70_df = pd.read_table(seq_id_70_clusters_filepath, sep='\t', header=None,
                                        names=['cluster-representative', 'cluster-member'])
    benchmark_metadata_df = pd.read_csv(benchmark_metadata_filepath)
    df = add_clusters(pairs_df=df, col1='protein1', col2='protein2', clusters_df=clusters_seqid70_df)
    df.rename(columns={'protein1_cluster': 'protein1_clusters_seqid70',
                       'protein2_cluster': 'protein2_clusters_seqid70',
                       'pair_cluster': 'pair_cluster_seqid70'}, inplace=True)
    df_deduplicated1 = df[
        ~df['pair_cluster_seqid70'].isin(benchmark_metadata_df['pair_cluster'].to_list())]
    print(df_deduplicated1.shape)

    # Remove duplicate pairs within the AP-MS datasets (based on 95% seq id clusters):
    clusters_seqid95_df = pd.read_table(seq_id_95_clusters_filepath, sep='\t', header=None,
                                        names=['cluster-representative', 'cluster-member'])
    df_deduplicated1 = add_clusters(pairs_df=df_deduplicated1, col1='protein1', col2='protein2',
                                    clusters_df=clusters_seqid95_df)

    df_deduplicated1.rename(columns={'protein1_cluster': 'protein1_cluster_seqid95',
                                     'protein2_cluster': 'protein2_cluster_seqid95',
                                     'pair_cluster': 'pair_cluster_seqid95'}, inplace=True)
    df_deduplicated2 = df_deduplicated1.drop_duplicates(subset=['pair_cluster_seqid95'], keep='first')
    print(df_deduplicated2.shape)

    df_deduplicated2.to_csv(output_filepath, index=False)


def remove_duplicates_apms_results(results_dir, deduplicated_metadata_filepath, output_filepath):
    results_files = glob.glob(os.path.join(results_dir, '*.csv'))

    dataframes = []
    exclude = ['EV71-HEK293T', 'IAV-HTBE']
    for res_file in results_files:
        if 'sars_mers_sars2_human_ppi' not in res_file:
            dataset = os.path.split(res_file)[-1].split('_')[0]
        else:
            dataset = 'sars_mers_sars2_human_ppi'

        if dataset not in exclude:
            df = pd.read_csv(res_file)
            df['dataset'] = dataset
            dataframes.append(df)
    concat_df = pd.concat(dataframes, axis=0)

    # Remove results for models with nan coordinates (empty AF2 scores columns):
    concat_df = concat_df.dropna(subset=['ptm', 'iptm', 'ranking_confidence'], axis=0)

    # Drop duplicate protein pairs (protein pairs that were tested in more than 1 cell line):
    concat_df = concat_df.sort_values('iptm', ascending=False).drop_duplicates(subset='name', keep='first')

    concat_df[['protein1', 'protein2']] = concat_df['name'].str.split(pat='_', expand=True)

    print(concat_df.shape)
    metadata_df = pd.read_csv(deduplicated_metadata_filepath)
    results_df_deduplicated = copy.deepcopy(concat_df[concat_df['name'].isin(metadata_df['protein_pair'].to_list())])
    print(results_df_deduplicated.shape)
    results_df_deduplicated.to_csv(output_filepath, index=False)


def remove_duplicates_apms_results_af3(results_filepath, deduplicated_metadata_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    print(results_df.shape)
    metadata_df = pd.read_csv(deduplicated_metadata_filepath)
    metadata_df['protein_pair'] = metadata_df['protein_pair'].str.upper().str.replace('^',
                                                                                      '-')  # because pairs had to be renamed when running AF3
    results_df_deduplicated = copy.deepcopy(results_df[results_df['name'].isin(metadata_df[
                                                                                   'protein_pair'].to_list())])  # unpublished pairs had already been removed from metadata_df, so this will remove them from the AF3 results file
    print(results_df_deduplicated.shape)
    results_df_deduplicated.to_csv(output_filepath, index=False)


if __name__ == '__main__':
    create_individual_fastas_hpidb(
        hpidb_dataset_path='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered.csv',
        output_dir='../data/hpidb2/fastafiles/individual_proteins')
    create_fasta_for_mmseqs(
        virus_benchmark_seqs_path='../data/virus_mammalia_dimers/fullseq_fasta_files/individual_proteins',
        bacteria_benchmark_seqs_path='../data/bacteria_mammalia_dimers/fullseq_fasta_files/individual_proteins',
        apms_seqs_path='../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
        sars_mers_sars2_path='../data/krogan_lab_host_pathogen_data/sars_mers_sars2_human_ppi/individual_proteins',
        hpidb_seqs_path='../data/hpidb2/fastafiles/individual_proteins',
        output_filepath='../data/all_seqs_min10aa_filenames_as_ids.fasta',
        minimum_length=10)
    apms_input_pairs('../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs.csv')

    add_benchmark_clusters(benchmark_pair_seqs_dir='../data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs',
                           seq_clusters_filepath='../data/clusterRes70_cluster.tsv',
                           output_filepath='../data/virus_mammalia_dimers_clusterRes70_groups.csv')
    add_benchmark_clusters(benchmark_pair_seqs_dir='../data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs',
                           seq_clusters_filepath='../data/clusterRes70_cluster.tsv',
                           output_filepath='../data/bacteria_mammalia_dimers_clusterRes70_groups.csv')

    remove_duplicates_fullseq_benchmark_metadata(
        bacteria_metadata_filepath='../data/bacteria_mammalia_dimers/70identity_search_results.csv',
        bacteria_clusters_filepath='../data/bacteria_mammalia_dimers_clusterRes70_groups.csv',
        virus_metadata_filepath='../data/virus_mammalia_dimers/70identity_search_results.csv',
        virus_clusters_filepath='../data/virus_mammalia_dimers_clusterRes70_groups.csv',
        output_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv')
    remove_duplicates_benchmark_results(
        results_filepath='../results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv',
        deduplicated_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        output_filepath='../results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_deduplicated.csv')
    remove_duplicates_benchmark_results(
        results_filepath='../results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv',
        deduplicated_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        output_filepath='../results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_deduplicated.csv')
    remove_duplicates_benchmark_results(
        results_filepath='../results/af3_benchmark/bacteria_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues.csv',
        deduplicated_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        output_filepath='../results/af3_benchmark/bacteria_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues_deduplicated.csv')
    remove_duplicates_benchmark_results(
        results_filepath='../results/af3_benchmark/virus_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues.csv',
        deduplicated_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        output_filepath='../results/af3_benchmark/virus_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues_deduplicated.csv')

    remove_duplicates_apms_input_pairs(
        input_pairs_filepath='../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs.csv',
        seq_id_70_clusters_filepath='../data/clusterRes70_cluster.tsv',
        benchmark_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        seq_id_95_clusters_filepath='../data/clusterRes95cov_cluster.tsv',
        output_filepath='../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs_deduplicated.csv')
    remove_duplicates_apms_results(results_dir='../results/krogan_lab_datasets/only_best_model',
                                   deduplicated_metadata_filepath='../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs_deduplicated.csv',
                                   output_filepath='../results/krogan_lab_datasets/only_best_model/all/results_published_datasets_onlybestmodel_deduplicated.csv')
    remove_duplicates_apms_results_af3(
        results_filepath="../results/krogan_lab_datasets/af3/krogan_lab_datasets_af3_results.csv",
        deduplicated_metadata_filepath='../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs_deduplicated.csv',
        output_filepath="../results/krogan_lab_datasets/af3/krogan_lab_datasets_af3_results_deduplicated.csv")

    remove_duplicates_hpidb_metadata(
        hpidb_metadata_filepath='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered.csv',
        seq_id_70_clusters_filepath='../data/clusterRes70_cluster.tsv',
        benchmark_metadata_filepath='../data/fullseq_benchmark_dataset_deduplicated.csv',
        seq_id_95_clusters_filepath='../data/clusterRes95cov_cluster.tsv',
        apms_filepath='../data/krogan_lab_host_pathogen_data/krogan_lab_input_pairs_deduplicated.csv',
        output_filepath='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered_deduplicated.csv')
    remove_duplicates_hpidb_results(results_filepath='../results/hpidb/hpidb_af2_results.csv',
                                    deduplicated_metadata_filepath='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered_deduplicated.csv',
                                    output_filepath='../results/hpidb/hpidb_filtered_af2_results_deduplicated.csv')
    remove_duplicates_hpidb_results(results_filepath='../results/hpidb/af3/hpidb_af3_results.csv',
                                    deduplicated_metadata_filepath='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered_deduplicated.csv',
                                    output_filepath='../results/hpidb/af3/hpidb_af3_results_deduplicated.csv')
