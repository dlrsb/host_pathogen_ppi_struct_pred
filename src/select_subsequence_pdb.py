import copy
import os
import argparse

import pandas as pd
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


def get_residues_from_sifts(model_pdb_id, chain_id, uniprot_id):
    sifts_df = pd.read_csv('../data/sifts_pdb_chain_uniprot.csv', header=1)
    model_sifts_data = sifts_df[sifts_df['PDB'] == model_pdb_id]
    print(chain_id)
    chain_sifts = model_sifts_data[(model_sifts_data['SP_PRIMARY'] == uniprot_id) & (model_sifts_data['CHAIN'] == chain_id)]
    chain_sifts.reset_index(drop=True, inplace=True)
    print(chain_sifts)
    if chain_sifts.shape[0] > 1:
        start = 1
        stop = 1
        for i, row in chain_sifts.iterrows():
            if i == 0:
                start = row['SP_BEG']
            elif row['SP_BEG'] < start:
                start = row['SP_BEG']

            if row['SP_END'] > stop:
                stop = row['SP_END']
    else:
        start = chain_sifts['SP_BEG'].item()
        stop = chain_sifts['SP_END'].item()

    return start-1, stop


def get_residues_by_aligning(model_chain_seqrecord, native_seqrecords_dict):
    native_seqrecord = native_seqrecords_dict[model_chain_seqrecord.id]
    blosum62_alphabet = substitution_matrices.load('BLOSUM62').alphabet
    seq = model_chain_seqrecord.seq
    for char in seq:
        if char not in blosum62_alphabet:
            seq = seq.replace(char, 'X')
    model_chain_seqrecord.seq = seq
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(native_seqrecord.seq, model_chain_seqrecord.seq)[0]
    print(alignment.format())
    subsequence_residues = alignment.aligned[1].tolist()
    if len(subsequence_residues) == 1:
        start = subsequence_residues[0][0]
        stop = subsequence_residues[0][1]
    else:
        start = -1
        stop = -1
        for i, subseq_res in enumerate(subsequence_residues):
            if (start == -1 or subseq_res[0] < start) and (subseq_res[1] - subseq_res[0] > 2): # to avoid artifacts at the start of the native sequence that may align to random areas in the model sequence
                start = subseq_res[0]

            if subseq_res[1] > stop and (subseq_res[1] - subseq_res[0] > 2):  # to avoid artifacts at the end of the native sequence that may align to random areas in the model sequence
                stop = subseq_res[1]

    return start, stop


def create_selection(model_pdb_filepath, subsequence_residues):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    model_structure = parser.get_structure(os.path.splitext(os.path.split(model_pdb_filepath)[-1])[0], model_pdb_filepath)[0]
    new_model_structure = copy.deepcopy(model_structure)
    for chain in model_structure:
        #print(chain)
        residues_to_keep = list(range(subsequence_residues[chain.id][0]+1, subsequence_residues[chain.id][1]+1))
        for residue in chain:
            if residue.id[1] not in residues_to_keep:
                del new_model_structure[chain.id][residue.id]

    io.set_structure(new_model_structure)
    pdb_dir, filename = os.path.split(model_pdb_filepath)
    #new_filename = '{}_selected.pdb'.format(filename.split('.')[0])
    new_filename = '{}.selected'.format(filename)
    io.save(os.path.join(pdb_dir, new_filename))
    return os.path.join(pdb_dir, new_filename)


def check_selection(new_model_filepath, subsequences):
    new_model_seqrecords = list(SeqIO.parse(new_model_filepath, 'pdb-atom'))
    for chain in new_model_seqrecords:
        chain_id = chain.id.split(':')[-1]
        print(chain_id)
        print(chain.seq)
        print(subsequences[chain_id])
        if chain.seq != subsequences[chain_id]:
            print('Sequence is not the same as the one found when aligning!')


def select_residues_pdb(model_pdb_filepath, model_input_fasta, native_pdb_filepath):
    model_pdb_id = os.path.splitext(os.path.split(model_input_fasta)[-1])[0].split('-')[0].lower()

    model_seqrecords = list(SeqIO.parse(model_input_fasta, 'fasta')) # using FASTA instead of reading sequence from the PDB file directly because sequences that start with 'X' were causing problems
    native_seqrecords = list(SeqIO.parse(native_pdb_filepath, 'pdb-atom'))
    native_seqrecords_dict = {native_seqrecord.id.split(':')[-1]: native_seqrecord for native_seqrecord in native_seqrecords}
    model_residues_to_keep = {}
    subsequences = {}
    for model_seqrecord in model_seqrecords:
        print(model_seqrecord.id)
        uniprot_id = model_seqrecord.description.split(' ')[-1]
        print(uniprot_id)
        if model_pdb_id.upper() in ['8DLY', '8DLQ', '7TUQ', '6T36']:
            start, stop = get_residues_by_aligning(model_chain_seqrecord=model_seqrecord,
                                                   native_seqrecords_dict=native_seqrecords_dict)
            model_residues_to_keep[model_seqrecord.id] = [start, stop]
            subsequences[model_seqrecord.id] = model_seqrecord.seq[start:stop]
        else:
            try:
                start, stop = get_residues_from_sifts(model_pdb_id, model_seqrecord.id, uniprot_id)
                model_residues_to_keep[model_seqrecord.id] = [start, stop]
                subsequences[model_seqrecord.id] = model_seqrecord.seq[start:stop]
            except Exception as e:
                print(e)
                start, stop = get_residues_by_aligning(model_chain_seqrecord=model_seqrecord,
                                         native_seqrecords_dict=native_seqrecords_dict)
                #print(start)
                #(stop)
                model_residues_to_keep[model_seqrecord.id] = [start, stop]
                subsequences[model_seqrecord.id] = model_seqrecord.seq[start:stop]

    print(model_residues_to_keep)

    new_model_filepath = create_selection(model_pdb_filepath=model_pdb_filepath,
                                          subsequence_residues=model_residues_to_keep)

    check_selection(new_model_filepath=new_model_filepath,
                    subsequences=subsequences)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Select the residues from a model structure that match the residues in the native PDB file')
    arg_parser.add_argument(
        '-m',
        dest='model_pdb_filepath',
        default=None,
        help='Path to PDB file containing the AF-multimer predicted structure',
    )
    arg_parser.add_argument(
        '-f',
        dest='model_input_fasta',
        default=None,
        help='Path to the input FASTA file',
    )
    arg_parser.add_argument(
        '-n',
        dest='native_pdb_filepath',
        default=None,
        help='Path to the PDB file containing the native (experimentally-determined) structure',
    )
    args = arg_parser.parse_args()
    select_residues_pdb(**vars(args))