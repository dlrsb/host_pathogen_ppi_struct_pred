import os
import pickle
import numpy as np
from Bio import AlignIO


def count_hits_a3m(filepath):
    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>') and 'chain_' not in line:
                count += 1
    return count


def count_hits_sto(filepath):
    alignment = AlignIO.read(filepath, 'stockholm')
    return len(alignment) - 1 # subtracting 1 so that the query sequence is not included in the # of hits


def get_num_msa_hits(results_dir):
    results = {}
    with open(os.path.join(results_dir, 'features.pkl'), 'rb') as f:
        features = pickle.load(f)
    try:
        results['num_alignments'] = features['num_alignments'].item()
    except ValueError as e:
        results['num_alignments'] = str(features['num_alignments'].tolist())
    del features

    msas_dir = os.path.join(results_dir, 'msas')
    chains = [x for x in os.listdir(msas_dir) if os.path.isdir(os.path.join(msas_dir, x))]
    for chain in chains:
        chain_dir = os.path.join(msas_dir, chain)
        alignment_files = os.listdir(chain_dir)
        for alignment_file in alignment_files:
            name, extension = os.path.split(alignment_file)[-1].split('.')
            try:
                if extension == 'a3m':
                    results['num_{db}_{chain}'.format(db=name, chain=chain)] = count_hits_a3m(os.path.join(chain_dir, alignment_file))
                else:
                    results['num_{db}_{chain}'.format(db=name, chain=chain)] = count_hits_sto(os.path.join(chain_dir, alignment_file))
            except Exception as e:
                print(e)
                results['num_{db}_{chain}'.format(db=name, chain=chain)] = np.nan

    return results
