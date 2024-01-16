import os
import re
import copy

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from get_full_sequences_uniprot import get_seq_from_uniprot


def split_datasets(output_dir):
    df = pd.read_excel('C:/Users/dlrba/Dropbox/igc_postdoc/projects/virus_host_ppi_structure_prediction/host_pathogen_ppi_struct_pred/data/krogan_lab_host_pathogen_data/Krogan_PPI_multipathogen.xlsx',
                       sheet_name='Krogan_PPI_multipathogen_FinalI',
                       keep_default_na=False)
    no_prefix_needed = ['Dengue-HEK293T', 'Zika-HEK293T', 'HPV-HEK293T']
    for dataset in df['dataset'].unique().tolist():
        if 'SARS' not in dataset and 'MERS' not in dataset:
            dataset_df = df[df['dataset'] == dataset]
            dataset_df['Bait'] = dataset_df['Bait'].str.replace(' ', '-')
            dataset_df['Bait'] = dataset_df['Bait'].str.replace('_', '-')
            if dataset not in no_prefix_needed:
                dataset_df['protein_pair'] = dataset_df['pathogen'] + '-' + dataset_df['Bait'] + '_' + dataset_df['Prey']
            else:
                dataset_df['protein_pair'] = dataset_df['Bait'] + '_' + dataset_df['Prey']
            dataset_df.to_csv(os.path.join(output_dir, f'{dataset}.csv'), index=False)


def create_cvb3_dataset():
    dataset_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Enteroviruses/CVB3_AP-MS_Supplemental/CVB3_AP-MS_Table_4_Supplemental.xlsx',
                             sheet_name='APMS for CV-B3',
                             nrows=7742)
    high_confidence_pairs_df = dataset_df[(dataset_df['mist'] >= 0.75) & (dataset_df['organism'] == 'Homo sapiens (Human)')]
    high_confidence_pairs_df['bait'] = high_confidence_pairs_df['bait'].str.replace('CVB_', 'CVB3-')
    high_confidence_pairs_df['protein_pair'] = high_confidence_pairs_df['bait'] + '_' + high_confidence_pairs_df['prey']
    high_confidence_pairs_df.to_csv('../data/krogan_lab_host_pathogen_data/datasets/CVB3.csv', index=False)


def filter_ev71_dataset():
    dataset_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/EV71-HEK293T.csv')
    preys_exclude = ['EV71_2B', 'EV71_2BC', 'EV71_2C', 'CVB2_2BC', 'EV71_3AB', 'EV71_VP1']
    baits_exclude = ['Vector', 'Vif']
    dataset_df_filtered1 = dataset_df[~dataset_df['Prey'].isin(preys_exclude)]
    dataset_df_filtered2 = dataset_df_filtered1[~dataset_df_filtered1['Bait'].isin(baits_exclude)]
    dataset_df_filtered2.to_csv('../data/krogan_lab_host_pathogen_data/datasets/EV71-HEK293T_filtered_full.csv', index=False)


def remove_mg132_ev71():
    dataset_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/EV71-HEK293T_filtered_full.csv')
    dataset_df['protein_pair'] = dataset_df['protein_pair'].str.replace('-MG132', '')
    dataset_df.to_csv('../data/krogan_lab_host_pathogen_data/datasets/EV71-HEK293T_filtered_full_nomg132.csv', index=False)


def get_prey_fasta_file(prey_id, fasta_output_dir):
    prey_output_fasta_path = os.path.join(fasta_output_dir, '{}.fasta'.format(prey_id))
    if not os.path.exists(prey_output_fasta_path):
        prey_record = get_seq_from_uniprot(
            prey_id)  # getting from Uniprot as the file in the OneDrive folder says it's from 2016
        prey_record.id = prey_id
        prey_record.name = ''
        prey_record.description = ''
        SeqIO.write(prey_record, prey_output_fasta_path, 'fasta-2line')


def create_fasta_files_chlamydia():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/Chlamydia-HEK293T.csv')
    fasta_records = list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Bacterial/Uniprot.Ctrachomatis.2019.04.03.fasta',
                               format='fasta'))

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join('../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins','{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            if bait == 'Chlamydia-CT442-71-150':
                bait_id = 'srp'
                residues = [71, 150]
            elif bait.startswith('Chlamydia-CT'):
                bait_id = bait[10:12] + '_' + bait[12:15]
                residues = [int(x) for x in bait[16:].split('-') if x]  # get residues, if they exist
            elif bait == 'Chlamydia-IncA':
                bait_id = 'incA'
                residues = []
            else:
                bait_id = bait[10:20].replace('-', '_')
                residues = [int(x) for x in bait[21:].split('-')[0:2] if x and x != 'MS']

            for record in fasta_records:
                if (f'GN={bait_id}' in record.description) or (bait_id in record.id):
                    bait_record = record
                    bait_record.id = bait
                    bait_record.name = ''
                    bait_record.description = ''
                    if len(residues) != 0:
                        subseq = bait_record.seq[residues[0]-1:residues[1]]
                        bait_record.seq = subseq
                    SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_mtb():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/Mtb-HEK293T.csv')
    pairs_df['bait_lower'] = pairs_df['Bait'].str.lower()
    mtb_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Bacterial/TB (2018).xlsx')
    mtb_df['gene_name_lower'] = mtb_df['Gene Name'].str.strip().str.lower()
    merged_df = pairs_df.merge(mtb_df, how='left', left_on='bait_lower', right_on='gene_name_lower')
    merged_df['Bait_id'] = merged_df['pathogen'] + '-' + merged_df['Bait']

    for index, row in merged_df.iterrows():
        pair = row['protein_pair']
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_uniprot_id = row['Uniprot']
        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins', '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_record = get_seq_from_uniprot(bait_uniprot_id)
            bait_record.id = bait
            bait_record.name = ''
            bait_record.description = ''
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_dengue():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/Dengue-HEK293T.csv')
    seqs_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Flaviviruses/DENV_ZIKV (2018).xlsx', sheet_name='Table S1', header=1)
    dengue_df = seqs_df[seqs_df['Bait'].str.contains('DENV2 16681')]
    dengue_df['Bait_renamed'] = dengue_df['Bait'].str.replace(' ', '-')
    tag = 'LEGGGGWSHPQFEKGGGSGGGSGGGSWSHPQFEKGPV*'

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_seq = dengue_df.loc[dengue_df['Bait_renamed'] == bait, 'Protein Sequence'].item()
            if bait_seq.endswith(tag):
                bait_seq = bait_seq[0:-len(tag)]
            bait_record = SeqRecord(Seq(bait_seq), id=bait, name='', description='')
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_zika():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/Zika-HEK293T.csv')
    seqs_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Flaviviruses/DENV_ZIKV (2018).xlsx', sheet_name='Table S1', header=1)
    zika_df = seqs_df[seqs_df['Bait'].str.contains('ZIKVfp')]
    zika_df['Bait_renamed'] = zika_df['Bait'].str.replace(' ', '-')
    tag = 'LEGGGGWSHPQFEKGGGSGGGSGGGSWSHPQFEKGPV*'

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_seq = zika_df.loc[zika_df['Bait_renamed'] == bait, 'Protein Sequence'].item()
            if bait_seq.endswith(tag):
                bait_seq = bait_seq[0:-len(tag)]
            bait_record = SeqRecord(Seq(bait_seq), id=bait, name='', description='')
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_ebola():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/Ebola-HEK293T.csv')
    seqs_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Krogan_PPI_multipathogen_FinalInteractomes_v2 20231126.xlsx', sheet_name='EBOV', header=12, nrows=7)
    seqs_df.rename(mapper={'Unnamed: 0': 'bait', 'Unnamed: 3': 'protein_seq'}, axis=1, inplace=True)
    seqs_df['bait_renamed'] = seqs_df['bait'].str.replace(' *', '')

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_seq = seqs_df.loc[seqs_df['bait_renamed'] == bait.split('-')[-1], 'protein_seq'].item()
            bait_record = SeqRecord(Seq(bait_seq), id=bait, name='', description='')
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_hpv():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/HPV-HEK293T.csv')
    seqs_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Krogan_PPI_multipathogen_FinalInteractomes_v2 20231126.xlsx', sheet_name='HPV', header=2, nrows=9)
    seqs_df.rename(mapper={'Unnamed: 0': 'bait', 'Amino acid sequence from SnapGene files': 'protein_seq'}, axis=1, inplace=True)
    seqs_df['bait_renamed'] = seqs_df['bait'].str.replace('_', '-')

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_seq = seqs_df.loc[seqs_df['bait_renamed'] == bait, 'protein_seq'].item()
            bait_record = SeqRecord(Seq(bait_seq), id=bait, name='', description='')
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_hbv():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/HBV-HUH7.csv')
    seqs_df = pd.read_excel('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Krogan_PPI_multipathogen_FinalInteractomes_v2 20231126.xlsx', sheet_name='HBV', header=2, nrows=8)
    seqs_df.rename(mapper={'Unnamed: 0': 'bait', 'Unnamed: 3': 'protein_seq'}, axis=1, inplace=True)
    seqs_df['bait_renamed'] = seqs_df['bait'].replace('POL', 'Pol')

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_seq = seqs_df.loc[seqs_df['bait_renamed'] == bait.split('-')[-1], 'protein_seq'].item()
            bait_record = SeqRecord(Seq(bait_seq), id=bait, name='', description='')
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_hiv():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/HIV-HEK293T.csv')

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

    # Only getting prey sequences as the bait sequences were manually obtained from a supplementary PDF file


def create_fasta_files_hcv():
    pairs_df1 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/HCV-HEK293T.csv')
    pairs_df2 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/HCV-HUH7.csv')

    for pair in pairs_df1['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

    for pair in pairs_df2['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

    # Only getting prey sequences as the bait sequences were manually obtained from a supplementary PDF file


def create_fasta_files_wnv():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/WNV-HEK293T.csv')

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        # Only creating FASTA files for the prey proteins in this function because the sequences for the WNV proteins
        # were manually obtained from the Uniprot entry for the WNV genome polyprotein


def create_fasta_files_cvb3():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/CVB3.csv')
    cvb3_fasta_records = {record.id.replace('_', '-').split('|')[0]: record  for record in list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Enteroviruses/CVB3 uniprot accession numbers for proteins.txt', 'fasta'))}
    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_record = cvb3_fasta_records[bait]
            bait_record.id = bait
            bait_record.name = ''
            bait_record.description = ''
            SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_kshv():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/KSHV-HEK293T.csv')
    kshv_fasta_records = list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/KSHV (HHV8)/SwissProt.HHV8.fasta', 'fasta')) + list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/KSHV (HHV8)/SwissProt.Human.HHV8.2016.01.11.fasta', 'fasta'))
    kshv_fasta_records = [record for record in kshv_fasta_records if 'OS=Homo sapiens' not in record.description]
    kshv_fasta_records_dict = {}
    for record in kshv_fasta_records:
        key = re.search('GN=(.+?) ', record.description).group()[3:-1]
        kshv_fasta_records_dict[key] = record

    baits_to_genes = {'K10': 'vIRF-4', 'K10.5': 'vIRF-3', 'K9': 'vIRF-1', 'ORF25': 'MCP', 'ORF26': 'TRX2',
                      'ORF47': 'gL', 'ORF6': 'DBP', 'ORF60': 'RIR2', 'ORF62': 'TRX1', 'ORF67': 'NEC2',
                      'ORF69': 'NEC1', 'ORF70': '70'}

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            bait_key = '.'.join(bait.split('-')[1:])
            if bait_key in baits_to_genes:
                bait_key = baits_to_genes[bait_key]
            if bait_key in kshv_fasta_records_dict.keys():
                bait_record = kshv_fasta_records_dict[bait_key]
                bait_record.id = bait
                bait_record.name = ''
                bait_record.description = ''
                SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_fasta_files_ev71():
    pairs_df = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/EV71-HEK293T_filtered_full_nomg132.csv')
    ev71_fasta_records = {record.id.replace('_', '-').split('|')[0]: record for record in list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/Enteroviruses/EV71 uniprot accession numbers for proteins.txt', 'fasta'))}

    for pair in pairs_df['protein_pair'].to_list():
        bait, prey = pair.split('_')

        get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        if bait in ev71_fasta_records.keys():
            bait_record = ev71_fasta_records[bait]
            bait_record.id = bait
            bait_record.name = ''
            bait_record.description = ''
            bait_output_fasta_path = os.path.join(
                '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
                '{}.fasta'.format(bait))
            if not os.path.exists(bait_output_fasta_path):
                SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')
        elif bait == 'EV71-VP1-2A-C110S':
            bait_prot1_record = ev71_fasta_records['EV71-VP1']
            bait_prot2_record = ev71_fasta_records['EV71-2A']
            mutated_seq = bait_prot2_record.seq[:110-1] + 'S' + bait_prot2_record.seq[110:]
            bait_record = SeqRecord(Seq(bait_prot1_record.seq+mutated_seq), id=bait, name='', description='')
            bait_output_fasta_path = os.path.join(
                '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
                '{}.fasta'.format(bait))
            if not os.path.exists(bait_output_fasta_path):
                SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')
        else: # point mutations
            bait_prot_name = '-'.join(bait.split('-')[0:-1])
            original_residue = bait.split('-')[-1][0]
            mutated_residue = bait.split('-')[-1][-1]
            if bait == 'EV71-P3-C147A':
                position = 108+int(bait.split('-')[-1][1:-1])-1
            else:
                position = int(bait.split('-')[-1][1:-1])-1
            bait_record = copy.deepcopy(ev71_fasta_records[bait_prot_name])
            seq = copy.deepcopy(bait_record.seq)
            if seq[position] == original_residue:
                mutated_seq = seq[:position] + mutated_residue + seq[position+1:]
                new_bait_record = SeqRecord(Seq(mutated_seq), id=bait, name='', description='')
                bait_output_fasta_path = os.path.join(
                    '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
                    '{}.fasta'.format(bait))
                if not os.path.exists(bait_output_fasta_path):
                    SeqIO.write(new_bait_record, bait_output_fasta_path, 'fasta-2line')
            else:
                print(bait)
                print(position)
                print(seq[position])
                print(original_residue)
                print(mutated_residue)


def create_fasta_files_iav():
    pairs_df1 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-A549.csv')
    pairs_df2 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-HTBE.csv')
    pairs_df3 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-THP1.csv')
    pairs_df4 = pd.read_csv('../data/krogan_lab_host_pathogen_data/datasets/IAV-HEK293T.csv')
    pairs_df_concat = pd.concat([pairs_df1, pairs_df2, pairs_df3, pairs_df4], axis=0, ignore_index=True)
    iav_fasta_records = {'IAV-'+record.id.replace('_', '-').split('|')[1]: record  for record in list(SeqIO.parse('../data/krogan_lab_host_pathogen_data/viral_bait_amino_acid_sequences/IAV_H1N1_H3N2_H5N1_ProteinSequences.fasta', 'fasta'))}

    for pair in pairs_df_concat['protein_pair'].to_list():
        bait, prey = pair.split('_')
        if pair != 'IAV-PB1F2_NA' and prey != 'eGFP':
            get_prey_fasta_file(prey, '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins')

        bait_output_fasta_path = os.path.join(
            '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
            '{}.fasta'.format(bait))
        if not os.path.exists(bait_output_fasta_path):
            if bait == 'IAV-eGFP-MOCK':
                bait = 'IAV-eGFP'
            if bait in iav_fasta_records.keys():
                bait_record = iav_fasta_records[bait]
                bait_record.id = bait
                bait_record.name = ''
                bait_record.description = ''
                SeqIO.write(bait_record, bait_output_fasta_path, 'fasta-2line')


def create_pair_fastafiles(dataset_dir, individual_fasta_dir, pair_fasta_dir):
    datasets = os.listdir(dataset_dir)
    for dataset in datasets:
        if dataset != 'EV71-HEK293T.csv':
            df = pd.read_csv(os.path.join(dataset_dir, dataset))
            for pair in df['protein_pair'].to_list():
                print('pair: '+ pair)
                bait, prey = pair.split('_')
                pair_fasta_path = os.path.join(pair_fasta_dir, f'{pair}.fasta')
                bait_path = os.path.join(individual_fasta_dir, f'{bait}.fasta')
                prey_path = os.path.join(individual_fasta_dir, f'{prey}.fasta')
                if os.path.exists(bait_path) and os.path.exists(prey_path) and not os.path.exists(pair_fasta_path):
                    bait_record = SeqIO.read(bait_path, 'fasta')
                    prey_record = SeqIO.read(prey_path, 'fasta')
                    SeqIO.write([bait_record, prey_record], os.path.join(pair_fasta_dir, f'{pair}.fasta'), 'fasta-2line')


if __name__ == '__main__':
    split_datasets('../../host_pathogen_ppi_struct_pred/data/krogan_lab_host_pathogen_data/datasets')
    create_cvb3_dataset()
    filter_ev71_dataset()
    remove_mg132_ev71()
    create_fasta_files_chlamydia()
    create_fasta_files_mtb()
    create_fasta_files_dengue()
    create_fasta_files_zika()
    create_fasta_files_ebola()
    create_fasta_files_hpv()
    create_fasta_files_hbv()
    create_fasta_files_hiv()
    create_fasta_files_hcv()
    create_fasta_files_wnv()
    create_fasta_files_cvb3()
    create_fasta_files_kshv()
    create_fasta_files_ev71()
    create_fasta_files_iav()
    create_pair_fastafiles('../data/krogan_lab_host_pathogen_data/datasets',
                           '../data/krogan_lab_host_pathogen_data/fastafiles/individual_proteins',
                           '../data/krogan_lab_host_pathogen_data/fastafiles/protein_pairs')
