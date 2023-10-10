import os
import subprocess
import pickle
from collections import defaultdict
import numpy as np
from Bio.PDB import Superimposer, CEAligner, PDBParser, PDBIO
from pymol import cmd
from pdockq_folddock import read_pdb, calc_pdockq
from pdockq_ifresid_huintaf2 import pDockQ, sigmoid
from pdockq2 import retrieve_IFplddt, retrieve_IFPAEinter, calc_pmidockq
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


def pdockq2_score(predicted_structure_filepath, results_pickle_filepath, cutoff=8):
    pdbp = PDBParser(QUIET=True)
    iopdb = PDBIO()

    structure = pdbp.get_structure('', predicted_structure_filepath)
    chains = []
    for chain in structure[0]:
        chains.append(chain.id)

    remain_contact_lst = []
    # retrieve interface plDDT at chain-level
    plddt_lst = []
    for idx in range(len(chains)):
        chain2_lst = list(set(chains) - set(chains[idx]))
        IF_plddt, contact_lst = retrieve_IFplddt(structure, chains[idx], chain2_lst, cutoff)
        plddt_lst.append(IF_plddt)
        remain_contact_lst.append(contact_lst)

    # retrieve interface PAE at chain-level
    with open(results_pickle_filepath, 'rb') as f:
        data = pickle.load(f)

    avgif_pae = retrieve_IFPAEinter(structure, data['predicted_aligned_error'], remain_contact_lst, cutoff)

    # calculate pmiDockQ
    res = calc_pmidockq(avgif_pae, plddt_lst)
    # res.index = chains
    res_dict = res.to_dict()
    res_dict_new = {}
    for k, score_dict in res_dict.items():
        for ch, v in score_dict.items():
            if k != 'prot':
                if k == 'pmidockq':
                    new_key = 'pdockq2_chain{chain}_{cutoff}Acutoff'.format(chain=ch, cutoff=cutoff)
                else:
                    new_key = 'pdockq2_{metric}_chain{chain}_{cutoff}Acutoff'.format(metric=k, chain=ch, cutoff=cutoff)
                res_dict_new[new_key] = v

    return res_dict_new


def rmsd_score_pymol(native_structure_filepath, predicted_structure_filepath, method_name='align',
                     outlier_rejection_cycles=5, atom_subset='all'):
    # TODO: figure out why super is not giving the same values as the PyMol GUI
    # https://bioinformatics.stackexchange.com/questions/19608/pymol-alignment-script
    # outlier_rejection_cycles=0 for all-atom RMSD
    methods_map = {'align': cmd.align,
                   'super': cmd.super,
                   'cealign': cmd.cealign}

    cmd.load(native_structure_filepath, 'experimental')
    cmd.load(predicted_structure_filepath, 'predicted')
    method = methods_map[method_name]  # TODO add a try except statement

    if method_name == 'align' or method_name == 'super':
        if atom_subset == 'CA':
            rmsd = method('predicted////CA', 'experimental////CA', cycles=outlier_rejection_cycles)[0]
        elif atom_subset == 'backbone':
            rmsd = method('predicted & backbone', 'experimental & backbone', cycles=outlier_rejection_cycles)[0]
        elif atom_subset == 'all':
            rmsd = method('predicted', 'experimental', cycles=outlier_rejection_cycles)[0]
        else:
            raise ValueError("Unknown atom_subset. Please choose one of 'CA', 'backbone', or 'all'")
    else:
        rmsd = method('experimental', 'predicted')['RMSD']
    return rmsd


def rmsd_score_biopython_superimposer(native_structure_filepath, predicted_structure_filepath, atom_subset='all'):
    # TODO: find a better way to get atom lists when atom_subset='all'
    # TODO: find a way to ignore atoms that are only present in one of the structures and only keep those present in both
    # atom_subset='all' is not working for at least 1 example (7NFR).
    # subset can be "CA", "backbone" or "all"
    # https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_superposition/
    # https://gist.github.com/JoaoRodrigues/e3a4f2139d10888c679eb1657a4d7080

    subsets = ['CA', 'backbone', 'all']
    if atom_subset not in subsets:
        raise ValueError("Unknown atom_subset. Please choose one of 'CA', 'backbone', or 'all'")

    native_structure = load_structure(native_structure_filepath)
    predicted_structure = load_structure(predicted_structure_filepath)

    # Align sequences for each chain
    alignments = align_seqs(native_structure_filepath, predicted_structure_filepath)

    chain_ids = alignments.keys()

    native_structure_atoms = []
    predicted_structure_atoms = []
    for chain_id in chain_ids:
        native_chain = native_structure[chain_id]
        predicted_chain = predicted_structure[chain_id]
        native_seq_aligned, predicted_seq_aligned, _, _, _ = alignments[chain_id]
        mask = [True if aa != '-' else False for aa in native_seq_aligned]

        native_all_atoms_ids = defaultdict(list)
        for i, residue in enumerate(native_chain):
            if atom_subset == 'CA':
                native_structure_atoms.append(residue['CA'])
            elif atom_subset == 'backbone':
                for atom_type in ['CA', 'C', 'N', 'O']:
                    native_structure_atoms.append(residue[atom_type])
            else:  # subset='all'
                for atom in residue:
                    native_structure_atoms.append(atom)
                    native_all_atoms_ids[str(i)].append((chain_id, residue.get_resname(), atom.id))

        j = -1
        for residue, use in zip(predicted_chain, mask):
            if use:
                j += 1
                if atom_subset == 'CA':
                    predicted_structure_atoms.append(residue['CA'])
                elif atom_subset == 'backbone':
                    for atom_type in ['CA', 'C', 'N', 'O']:
                        predicted_structure_atoms.append(residue[atom_type])
                else:  # subset='all'
                    for atom in residue:
                        if (chain_id, residue.get_resname(), atom.id) in native_all_atoms_ids[str(j)]:
                            # only include the atom if it's also part of the same residue in the ground truth structure
                            predicted_structure_atoms.append(atom)
    sup = Superimposer()
    sup.set_atoms(native_structure_atoms, predicted_structure_atoms)  # sup.set_atoms(fixed, moving)
    rmsd = sup.rms

    return rmsd


def rmsd_score_biopython_cealign(native_structure, predicted_structure, window_size=8, max_gap=30):
    cealigner = CEAligner(window_size=window_size, max_gap=max_gap)
    cealigner.set_reference(native_structure)
    cealigner.align(predicted_structure, transform=False)
    rmsd = cealigner.rms
    return rmsd


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


def dockq_score(native_structure_filepath, predicted_structure_filepath):
    predicted_structure = load_structure(predicted_structure_filepath)
    chains = [chain.id for chain in predicted_structure]
    native_structure = load_structure(native_structure_filepath)
    native_chains_in_predicted_struct = [chain.id for chain in native_structure if chain.id in chains]
    # if chains != native_chains_in_predicted_struct: # chains have different order in the native_structure file
    #     # sort the native_structure_file so that DockQ doesn't warn about a chain mismatch
    #     try:
    #         split_filepath ='{}_1{}'.format(os.path.splitext(native_structure_filepath)[0],
    #                                         os.path.splitext(native_structure_filepath)[1])
    #         sorted_filepath = '{}_sorted{}'.format(os.path.splitext(native_structure_filepath)[0],
    #                                         os.path.splitext(native_structure_filepath)[1])
    #         wdir = os.path.split(native_structure_filepath)[0]
    #         print(wdir)
    #         # TODO: create a small bash script to run this
    #         subprocess.run(
    #         ['conda', 'activate', 'vh_struct_pred', '&&', 'cd', wdir, '&&', 'pdb_splitmodel', '-C',
    #          native_structure_filepath,
    #          '&&', 'pdb_sort', split_filepath, '>',
    #          sorted_filepath, '&&', 'rm', split_filepath], shell=True, check=True)
    #         native_structure_filepath = sorted_filepath
    #     except subprocess.CalledProcessError as err:
    #         print(err)

    if len(chains) == 2:
        command = ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_dockq.sh', '-m',
                   predicted_structure_filepath, '-n', native_structure_filepath, '--native_chain1', chains[0],
                   '--model_chain1', chains[0], '--native_chain2', chains[1], '--model_chain2', chains[1]]
    else:
        # TODO: need to add the commands for protein complexes with more than 2 chains (calculate DockQ for all pairs
        #  of chains
        pass

    if chains != native_chains_in_predicted_struct:
        command += ['--sort_native']

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
    scores_to_keep = ["ptm", "iptm", "ranking_confidence"]
    with open(results_pickle_filepath, 'rb') as f:
        af_results = pickle.load(f)
    selected_scores = {k: af_results[k].item() for k in scores_to_keep}
    return selected_scores
