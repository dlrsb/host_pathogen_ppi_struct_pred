import os
import json
import pandas as pd
from Bio import SeqIO


def create_json(protein_pairs_dir, output_dir, seed, use_templates=False):
    protein_pairs_filepaths = [os.path.join(protein_pairs_dir, x) for x in os.listdir(protein_pairs_dir)]
    for pair_fasta in protein_pairs_filepaths:
        pair_name = os.path.splitext(os.path.split(pair_fasta)[-1])[0].lower()
        records = list(SeqIO.parse(pair_fasta, 'fasta'))
        chains = [record.id for record in records]

        seq1 = records[0].seq
        seq2 = records[1].seq
        if len(chains[0]) > 1:
            id1 = 'A'
            id2 = 'B'
        else:
            id1 = records[0].id
            id2 = records[1].id

        if use_templates:
            job_str = '{"name": "%s", "sequences": [{"protein": {"id": "%s", "sequence": "%s"}}, {"protein": {"id": "%s", "sequence": "%s"}}], "modelSeeds": [%s], "dialect": "alphafold3", "version": 1}' % (pair_name, id1, seq1, id2, seq2, seed)
        else:
            job_str = '{"name": "%s", "sequences": [{"protein": {"id": "%s", "sequence": "%s", "templates": []}}, {"protein": {"id": "%s", "sequence": "%s", "templates": []}}], "modelSeeds": [%s], "dialect": "alphafold3", "version": 1}' % (pair_name, id1, seq1, id2, seq2, seed)
        json_object = json.loads(job_str)
        json_output_path = os.path.join(output_dir, f'{pair_name}.json')
        with open(json_output_path, 'w') as f:
            json.dump(json_object, f, indent=4)


def create_json_hpidb(protein_pairs_filepath, output_dir, seed, use_templates=False):
    df = pd.read_csv(protein_pairs_filepath)
    df['id'] = df['id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace('NP_', 'NP-')

    for i, r in df.iterrows():
        pair_name = r['id'].lower()
        split_seqs = r['sequence'].split(':')
        seq1 = split_seqs[0]
        seq2 = split_seqs[1]
        id1 = 'A'
        id2 = 'B'

        if use_templates:
            job_str = '{"name": "%s", "sequences": [{"protein": {"id": "%s", "sequence": "%s"}}, {"protein": {"id": "%s", "sequence": "%s"}}], "modelSeeds": [%s], "dialect": "alphafold3", "version": 1}' % (pair_name, id1, seq1, id2, seq2, seed)
        else:
            job_str = '{"name": "%s", "sequences": [{"protein": {"id": "%s", "sequence": "%s", "templates": []}}, {"protein": {"id": "%s", "sequence": "%s", "templates": []}}], "modelSeeds": [%s], "dialect": "alphafold3", "version": 1}' % (pair_name, id1, seq1, id2, seq2, seed)

        json_object = json.loads(job_str)
        json_output_path = os.path.join(output_dir, f'{pair_name}.json')
        with open(json_output_path, 'w') as f:
            json.dump(json_object, f, indent=4)


if __name__ == '__main__':
    create_json('../data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs', '../data/virus_mammalia_dimers/af3_input', '4')
    create_json('../data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs', '../data/bacteria_mammalia_dimers/af3_input', '4')
    create_json('../data/krogan_lab_host_pathogen_data/fastafiles/protein_pairs',
                '../data/krogan_lab_host_pathogen_data/af3_input', seed='4', use_templates=True)
    create_json('../data/krogan_lab_host_pathogen_data/sars_mers_sars2_human_ppi/protein_pairs',
                '../data/krogan_lab_host_pathogen_data/af3_input', seed='4', use_templates=True)
    create_json_hpidb('../data/hpidb2/hpidb_filtered/hpidb_selected_pairs_seqs_filtered.csv',
                '../data/hpidb2/hpidb_filtered/af3_input', seed='4', use_templates=True)

