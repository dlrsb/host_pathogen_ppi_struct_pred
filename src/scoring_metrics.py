import os
import subprocess
import pickle
from collections import defaultdict
import numpy as np
from Bio.PDB import PDBParser
from pdockq_folddock import read_pdb, calc_pdockq
from pdockq_ifresid_huintaf2 import pDockQ, sigmoid
from utils import load_structure, align_seqs


def pdockq_folddock_score(predicted_structure_filepath, cutoff=5):
    results = {}
    chain_coords, chain_plddt = read_pdb(predicted_structure_filepath)
    results['folddock_pdockq_{}Acutoff'.format(cutoff)], results['folddock_ppv_{}Acutoff'.format(cutoff)], results[
        'folddock_num_interface_contacts_{}Acutoff'.format(cutoff)] = calc_pdockq(chain_coords, chain_plddt,
                                                                                  cutoff)
    return results


def pdockq_huintaf2_score(predicted_structure_filepath, cutoff=5, factor=1.0):
    results = {}
    if cutoff <= 5:
        popt = [6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
    elif cutoff <= 6:
        popt = [7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
    elif cutoff <= 7:
        popt = [7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
    elif cutoff <= 8:
        popt = [7.18442739e-01, 3.60791204e+02, 3.01635944e-02, 2.04076969e-02]
    elif cutoff <= 9:
        popt = [7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
    elif cutoff <= 10:
        popt = [7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]  # used in Interaction studies.
    elif cutoff <= 11:
        popt = [7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
    else:
        popt = [7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]

    parser = PDBParser()
    structure = parser.get_structure(os.path.split(predicted_structure_filepath)[-1].split('.')[0],
                                     predicted_structure_filepath)
    NumRes, IF_plDDT, plDDT1, plDDT2, Diso1, Diso2, len1, len2, residues_seq1, residues_seq2 = pDockQ(
        structure, cutoff)
    tiny = 1.e-20
    results['huintaf2_pdockq_{}Acutoff'.format(cutoff)] = sigmoid(np.log(NumRes + tiny) * IF_plDDT * factor, *popt)
    results['huintaf2_num_residues_{}Acutoff'.format(cutoff)] = NumRes
    results['huintaf2_interface_plddt_{}Acutoff'.format(cutoff)] = IF_plDDT
    results['huintaf2_interface_residues_chain1_{}Acutoff'.format(cutoff)] = str(residues_seq1)
    results['huintaf2_interface_residues_chain2_{}Acutoff'.format(cutoff)] = str(residues_seq2)
    return results


def tm_score(native_structure_filepath, predicted_structure_filepath):
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/tools/MMalign', predicted_structure_filepath,
             native_structure_filepath, '-outfmt',
             '2'], check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        results_all = dict(zip(header, values))

        results = {}
        for new_key, old_key in zip(['mmalign_tmscore', 'mmalign_rmsd'], ['TM2', 'RMSD']):
            results[new_key] = results_all[old_key]
    except subprocess.CalledProcessError as err:
        print(err)
        results = {'mmalign_tmscore': np.nan, 'mmalign_rmsd': np.nan}

    return results


def dockq_score(native_structure_filepath, predicted_structure_filepath, fasta_filepath, select_residues=False):
    predicted_structure = load_structure(predicted_structure_filepath)
    chains = [chain.id for chain in predicted_structure]
    print(chains)
    native_structure = load_structure(native_structure_filepath)
    native_chains_in_predicted_struct = [chain.id for chain in native_structure if chain.id in chains]
    print(native_chains_in_predicted_struct)

    if len(chains) == 2:
        command = ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_dockq.sh', '-m',
                   predicted_structure_filepath, '-f', fasta_filepath, '-n', native_structure_filepath, '--native_chain1', chains[0],
                   '--model_chain1', chains[0], '--native_chain2', chains[1], '--model_chain2', chains[1]]
    else:
        # TODO: need to add the commands for protein complexes with more than 2 chains (calculate DockQ for all pairs
        #  of chains
        pass

    if chains != native_chains_in_predicted_struct:
        command += ['--sort_chains']

    if select_residues:
        command += ['--select_residues']

    try:
        called_process = subprocess.run(command, check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        dockq_line = None
        dockq_columns = ['DockQ', 'Fnat', 'iRMS', 'LRMS', 'Fnonnat']
        for line in lines:
            print(line)
            if all([x in line for x in dockq_columns]):
                dockq_line = line
                break
        if dockq_line is None:
            raise Exception('DockQ score not found in output!')
    except subprocess.CalledProcessError as err:
        print(err)
        results = {'dockq': np.nan, 'dockq_fnat': np.nan, 'dockq_irms': np.nan, 'dockq_lrms': np.nan,
                   'dockq_fnonnat': np.nan}
    except Exception as e:
        print(e)
        results = {'dockq': np.nan, 'dockq_fnat': np.nan, 'dockq_irms': np.nan, 'dockq_lrms': np.nan,
                   'dockq_fnonnat': np.nan}
    else:
        dockq_results = dockq_line.split(' ')
        keys = [k.lower() if k == 'DockQ' else 'dockq_{}'.format(k.lower()) for k in dockq_results[0:9:2]]
        values = dockq_results[1:10:2]
        results = dict(zip(keys, values))

    return results


def parse_dockq(result_filepath):
    with open(result_filepath, 'r') as f:
        result = [line.strip() for line in f.readlines()][-1].split(' ')[0:-2]
    scores_dict = {}
    for i in range(0, len(result), 2):
        scores_dict[result[i]] = float(result[i + 1])

    return scores_dict


def parse_mmalign(results_filepath):
    with open(results_filepath, 'r') as f:
        result = [line.strip() for line in f.readlines()]
    values = result[1].split('\t')
    scores_dict = {'mmalign_tmscore': values[3], 'mmalign_RMSD': values[4]}
    return scores_dict


def get_alphafold_results_scores(results_pickle_filepath):
    # read results pickle file and extract scores
    scores_to_keep = ['ptm', 'iptm', 'ranking_confidence']
    with open(results_pickle_filepath, 'rb') as f:
        af_results = pickle.load(f)
    selected_scores = {k: af_results[k].item() for k in scores_to_keep}
    return selected_scores
