import sys
import os
import subprocess
import copy
import requests
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, Select, PDBIO, Selection, NeighborSearch
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response


def load_structure(pdb_filepath):
    parser = PDBParser(QUIET=True)
    name = os.path.split(pdb_filepath)[-1].split('.')[0]
    structure = parser.get_structure(name, pdb_filepath)[0] # gets first model
    return structure


def run_pymol(cmds):
    #s_ = subprocess.run('module load gcc/6.3.0 pymol; pymol -cpQ', input='\n'.join(cmds), text=True, shell=True, capture_output=True)
    s_ = subprocess.run('pymol -cpQ', input='\n'.join(cmds), text=True, shell=True, capture_output=True)
    print(s_)
    return s_


def align_seqs(structure1_filepath, structure2_filepath):
    # based on https://gist.github.com/JoaoRodrigues/e3a4f2139d10888c679eb1657a4d7080
    # get sequences
    struct1_seqs = SeqIO.to_dict(SeqIO.parse(structure1_filepath, 'pdb-atom'))
    struct2_seqs = SeqIO.to_dict(SeqIO.parse(structure2_filepath, 'pdb-atom'))

    # align chains based on id (need to rename chains in 1 of the files beforehand)
    chains = struct1_seqs.keys()
    alignments = {}
    for chain in chains:
        seq1 = struct1_seqs[chain].seq.replace('X', '')
        seq2 = struct2_seqs[chain].seq.replace('X', '')
        alignment = pairwise2.align.globalds(
            seq1,
            seq2,
            substitution_matrices.load('BLOSUM62'),
            one_alignment_only=True,
            open=-10.0,
            extend=-0.5,
            penalize_end_gaps=(False, False),
        )
        alignments[chain.split(':')[-1]] = alignment[0]

    return alignments


def run_tmalign(model_structure_filepath, reference_structure_filepath):
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/tools/TMalign', model_structure_filepath,
             reference_structure_filepath, '-outfmt',
             '2'], check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        print(lines)
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        results_all = dict(zip(header, values))
        print(results_all)
        tmscore = results_all['TM2'] # using TM2, which is TM-score normalized by the reference structure
        rmsd = results_all['RMSD']
    except subprocess.CalledProcessError as err:
        print(err)
        tmscore = np.nan
        rmsd = np.nan
    except KeyError as err2:
        print(err2)
        tmscore = np.nan
        rmsd = np.nan

    return tmscore, rmsd


class SelectChains(Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_letters


class SelectResidues(Select):
    def __init__(self, chain_id, residues_list):
        self.chain_id = chain_id
        self.residues_list = residues_list

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        return residue.get_id()[1] in self.residues_list


def create_temp_pdb(pdb_filepath, chain_id, output_dir):
    struct = load_structure(pdb_filepath)
    io = PDBIO()
    pdb_name = os.path.split(os.path.splitext(pdb_filepath)[0])[-1]
    output_path = os.path.join(output_dir, f'{pdb_name}_chain{chain_id}.pdb')
    io.set_structure(struct)
    io.save(output_path, select=SelectChains(chain_id))
    return output_path


def create_interface_temp_pdb(pdb_filepath, chain_id, selected_residues, output_dir):
    struct = load_structure(pdb_filepath)
    io = PDBIO()
    pdb_name = os.path.split(os.path.splitext(pdb_filepath)[0])[-1]
    output_path = os.path.join(output_dir, f'{pdb_name}_chain{chain_id}_interface_residues.pdb')
    io.set_structure(struct)
    io.save(output_path, select=SelectResidues(chain_id, selected_residues))
    return output_path


def select_residues_within_radius(pdb_filepath, chain_id, residue_list, radius):
    struct = load_structure(pdb_filepath)
    atoms = Selection.unfold_entities(struct, 'A')
    ns = NeighborSearch(atoms)
    selected_residues = copy.deepcopy(residue_list)
    for res in residue_list:
        target_atom = struct[chain_id][res]['CA']
        close_residues = ns.search(target_atom.coord, radius, level='R')
        selected_residues.extend([x.get_id()[1] for x in close_residues])
    selected_residues_final = sorted(list(set(selected_residues)))
    return selected_residues_final


def calc_avg_plddt(pdb_filepath, interactor_chain, selected_residues):
    print(pdb_filepath)
    print(selected_residues)
    struct = load_structure(pdb_filepath)
    residue_plddt_vals = []
    for chain in struct:
        if chain.id == interactor_chain:
            interactor_residues = Selection.unfold_entities(chain, 'R')
            for residue in interactor_residues:
                if residue.id[1] in selected_residues:
                    b_factors_all_atoms = [x.get_bfactor() for x in Selection.unfold_entities(residue, 'A')]
                    print(b_factors_all_atoms[0])
                    residue_plddt_vals.append(Selection.unfold_entities(residue, 'A')[0].get_bfactor()) # B factor = plddt per residue. The value is the same for every atom in a residue

    # calculate mean plddt
    if len(residue_plddt_vals) == 0:
        avg_plddt = np.nan
        print('failed: '+pdb_filepath)
    else:
        avg_plddt = sum(residue_plddt_vals) / len(residue_plddt_vals)

    return avg_plddt


def get_interface_secondary_struct(pdb_filepath, chain, interface_residues):
    dssp_tuple = dssp_dict_from_pdb_file(pdb_filepath, DSSP='mkdssp', dssp_version='3.0.0')
    dssp_dict = dssp_tuple[0]
    if_ss_dict = {str(k[1][1]): v[1] for k, v in dssp_dict.items() if k[0] == chain and k[1][1] in interface_residues}
    counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
    for key, value in if_ss_dict.items():
        counts[value] += 1
    return counts


def sort_id(x):
    split_id = x.split('_')
    sorted_id = '_'.join(sorted(split_id))
    return sorted_id


def create_pymol_session(pdb_filepath1, shared_prot_chain1, interaction_id1, pdb_filepath2, shared_prot_chain2,
                         interaction_id2, output_filepath):
    cmd_ = ['delete all',
            'bg_color white',
            'space cmyk',
            f'load {pdb_filepath1}, {interaction_id1}',
            f'load {pdb_filepath2}, {interaction_id2}',
            f'color tv_orange, {interaction_id1}',
            f'color skyblue, {interaction_id2}',
            f'color gray50, {interaction_id1} & chain {shared_prot_chain1}',
            f'color gray50, {interaction_id2} & chain {shared_prot_chain2}',
            f'super {interaction_id1} & chain {shared_prot_chain1}, {interaction_id2} & chain {shared_prot_chain2}',
            f'save {output_filepath}']
    print(cmd_)
    run_pymol(cmd_)
