import os
import glob
import ast
import subprocess
from collections import Counter
import numpy as np
from scipy.stats import mode
import pandas as pd
from Bio.PDB import Selection
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

from utils import load_structure
from shared_interface_analysis import run_pymol


def calc_avg_plddt(pdb_filepath, interactor_chain, interface_residues):
    struct = load_structure(pdb_filepath)
    residue_plddt_vals = []
    for chain in struct:
        if chain.id == interactor_chain:
            interactor_residues = Selection.unfold_entities(chain, 'R')
            for residue in interactor_residues:
                if str(residue.id[1]) in interface_residues:
                    residue_plddt_vals.append(Selection.unfold_entities(residue, 'A')[0].get_bfactor()) # B factor = plddt per residue. The value is the same for every atom in a residue

    # calculate mean plddt
    avg_plddt = sum(residue_plddt_vals) / len(residue_plddt_vals)

    return avg_plddt


def create_pymol_session(pdb_filepath1, shared_prot_chain1, interaction_id1, pdb_filepath2, shared_prot_chain2,
                         interaction_id2, output_filepath):
    cmd_ = ['delete all',
            'bg_color white',
            'space cmyk',
            f'load {pdb_filepath1}, {interaction_id1}',
            f'load {pdb_filepath2}, {interaction_id2}',
            f'color gray70, {interaction_id1} & chain {shared_prot_chain1}',
            f'color gray70, {interaction_id2} & chain {shared_prot_chain2}',
            f'align {interaction_id1} & chain {shared_prot_chain1}, {interaction_id2} & chain {shared_prot_chain2}',
            f'save {output_filepath}']
    run_pymol(cmd_)


def run_ialign(structure1_filepath, structure2_filepath, scoring_metric, normalization_method, if_outputdir,
               distance_cutoff=4.5):
    called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_ialign2.sh', '-s1',
             structure1_filepath, '-s2', structure2_filepath, '-n', normalization_method, '-m', scoring_metric,
             '--dc', str(distance_cutoff), '--minp', '20', '--mini', '8', '--if_outputdir',
             if_outputdir],
            check=True)
    return called_process


def check_if_linear_old(interface_residues):
    diff = []
    for i in range(len(interface_residues)):
        if i != len(interface_residues) - 1:
            diff.append(interface_residues[i+1] - interface_residues[i])
    print(sum(diff) / len(diff))
    print(np.median(diff))
    counts = Counter(diff)
    print(counts[1]/len(diff))
    print(mode(diff))


def get_interface_secondary_struct(pdb_filepath, chain, interface_residues):
    dssp_tuple = dssp_dict_from_pdb_file(pdb_filepath, DSSP='mkdssp', dssp_version='3.0.0')
    dssp_dict = dssp_tuple[0]
    if_ss_dict = {str(k[1][1]): v[1] for k, v in dssp_dict.items() if k[0] == chain and k[1][1] in interface_residues}
    print(if_ss_dict)
    counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
    for key, value in if_ss_dict.items():
        counts[value] += 1
    print(counts)
    return counts


def analyze_top_scores(output_dir, only_top_human_interactor=False):
    print('filtering scores')
    df = pd.read_csv('../results/interface_analysis2/interface_comparison_results_all2.csv')
    #df_selected = df[(df['ialign8_itmscore'] >= 0.5) & (df['ialign8_isscore'] >= 0.3)]
    df_selected = df[(df['ialign8_itmscore'] >= 0.5) & (df['jaccard_index'] >= 0.3)]
    # df_selected.to_csv('../results/interface_analysis2/interface_comparison_results_selected.csv')
    if only_top_human_interactor:
        df_grouped_max = df_selected.loc[df_selected.groupby(by=['pathogen_interactor', 'bait'])['ialign8_itmscore'].idxmax()]
        df_selected = df_grouped_max
        #df_grouped_max.to_csv('../results/interface_analysis2/interface_comparison_results_selected_max.csv', index=False)
    metadata_df = pd.read_csv('../results/interface_analysis2/df_interactors_selected.csv')

    print('merging with metadata df')
    merged_df1 = df_selected.merge(metadata_df, how='left', left_on=['bait', 'pathogen_interactor'],
                                   right_on=['bait_id', 'interactor_id'])
    merged_df1.rename({'bait_id': 'pathogen_ppi_bait_id', 'bait_ifresid': 'pathogen_ppi_bait_ifresid',
                       'bait_chain': 'pathogen_ppi_bait_chain', 'interactor_chain': 'pathogen_ppi_interactor_chain',
                       'pdb': 'pathogen_ppi_pdb', 'interactor_ifresid': 'pathogen_ppi_interactor_ifresid',
                       'interaction_id': 'pathogen_interaction_id', 'interactor_id': 'pathogen_interactor_id'},
                      axis=1, inplace=True)
    merged_df1.drop(['pdockq', 'interaction_type'], axis=1, inplace=True)
    merged_df2 = merged_df1.merge(metadata_df, how='left', left_on=['bait', 'human_interactor'], right_on=['bait_id', 'interactor_id'])
    merged_df2.rename({'bait_id': 'human_ppi_bait_id', 'bait_ifresid': 'human_ppi_bait_ifresid',
                       'bait_chain': 'human_ppi_bait_chain', 'interactor_chain': 'human_ppi_interactor_chain',
                       'pdb': 'human_ppi_pdb', 'interactor_ifresid': 'human_ppi_interactor_ifresid',
                       'interaction_id': 'human_interaction_id', 'interactor_id': 'human_interactor_id'},
                      axis=1, inplace=True)
    merged_df2.drop(['pdockq', 'interaction_type'], axis=1, inplace=True)

    avg_plddt_dict = {'bait': [],
                      'pathogen_interactor': [],
                      'human_interactor': [],
                      'avg_if_plddt': [],
                      'full_struct_filepath': [],
                      'interface_filepath': []}

    ss_dict = {'bait': [],
               'pathogen_interactor': [],
               'human_interactor': [],
               'pathogen_H': [],
               'pathogen_B': [],
               'pathogen_E': [],
               'pathogen_G': [],
               'pathogen_I': [],
               'pathogen_T': [],
               'pathogen_S': [],
               'pathogen_-': [],
               'human_H': [],
               'human_B': [],
               'human_E': [],
               'human_G': [],
               'human_I': [],
               'human_T': [],
               'human_S': [],
               'human_-': [],
               'num_pathogen_if_residues': [],
               'num_human_if_residues': []}

    for index, row in merged_df2.iterrows():
        print(row['pathogen_interactor'])
        print(row['human_interactor'])
        # Separate cases by average interface plDDT (different output_dir)
        avg_plddt = calc_avg_plddt(pdb_filepath=row['pathogen_ppi_pdb'], interactor_chain=row['pathogen_ppi_interactor_chain'],
                       interface_residues=row['pathogen_ppi_interactor_ifresid'])
        print(avg_plddt)
        avg_plddt_dict['bait'].append(row['bait'])
        avg_plddt_dict['pathogen_interactor'].append(row['pathogen_interactor'])
        avg_plddt_dict['human_interactor'].append(row['human_interactor'])
        avg_plddt_dict['avg_if_plddt'].append(avg_plddt)

        if avg_plddt < 50.0:
            output_dir2 = os.path.join(output_dir, 'very_low_if_plddt')
        elif (avg_plddt >= 50.0) and (avg_plddt < 70.0):
            output_dir2 = os.path.join(output_dir, 'low_if_plddt')
        elif (avg_plddt >= 70.0) and (avg_plddt < 90.0):
            output_dir2 = os.path.join(output_dir, 'high_if_plddt')
        else:
            output_dir2 = os.path.join(output_dir, 'very_high_if_plddt')

        full_struct_output_filepath = os.path.join(output_dir2, f'full_{row["pathogen_interactor"]}_{row["human_interactor"]}_bait{row["pathogen_ppi_bait_id"]}.pse')
        interface_output_filepath = os.path.join(output_dir2, f'interface_{row["pathogen_interactor"]}_{row["human_interactor"]}_bait{row["pathogen_ppi_bait_id"]}.pse')
        avg_plddt_dict['full_struct_filepath'].append(full_struct_output_filepath)
        avg_plddt_dict['interface_filepath'].append(interface_output_filepath)

        # Determine secondary structure of pathogen interface:
        ss_dict['bait'].append(row['bait'])
        ss_dict['pathogen_interactor'].append(row['pathogen_interactor'])
        ss_dict['human_interactor'].append(row['human_interactor'])
        ss_dict['num_pathogen_if_residues'].append(len(ast.literal_eval(row['pathogen_ppi_interactor_ifresid'])))
        ss_dict['num_human_if_residues'].append(len(row['human_ppi_interactor_ifresid'].split(',')))
        pathogen_ss = get_interface_secondary_struct(pdb_filepath=row['pathogen_ppi_pdb'], chain=row['pathogen_ppi_interactor_chain'],
                                                     interface_residues=ast.literal_eval(row['pathogen_ppi_interactor_ifresid']))
        human_ss = get_interface_secondary_struct(pdb_filepath=row['human_ppi_pdb'], chain=row['human_ppi_interactor_chain'],
                                                  interface_residues=[int(x) for x in row['human_ppi_interactor_ifresid'].split(',')])
        for key in pathogen_ss:
            ss_dict[f'pathogen_{key}'].append(pathogen_ss[key])
        for key in human_ss:
            ss_dict[f'human_{key}'].append(human_ss[key])

        # Save PyMol session with full structures
        if not os.path.exists(full_struct_output_filepath):
            print('saving full structures')
            create_pymol_session(pdb_filepath1=row['pathogen_ppi_pdb'], shared_prot_chain1=row['pathogen_ppi_bait_chain'],
                                 interaction_id1=row['pathogen_interaction_id'], pdb_filepath2=row['human_ppi_pdb'],
                                 shared_prot_chain2=row['human_ppi_bait_chain'], interaction_id2=row['human_interaction_id'],
                                 output_filepath=full_struct_output_filepath)

        # Run ialign with selected metric and get interface PDB file (PDB files with _int, in "outputs" folder)
        if not os.path.exists(interface_output_filepath):
            print('running ialign to get interface PDB files')
            interface_output_dir = '/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/ialign_interface_pdbs'
            run_ialign(structure1_filepath=row['pathogen_ppi_pdb'], structure2_filepath=row['human_ppi_pdb'],
                       scoring_metric='tm', normalization_method='ave', distance_cutoff=8,
                       if_outputdir=interface_output_dir)
            interface_pdb_files = glob.glob(os.path.join(interface_output_dir, '*_int.pdb'))
            for pdb_file in interface_pdb_files:
                if row['human_interaction_id'] in pdb_file:
                    human_if_pdb_path = pdb_file
                else:
                    pathogen_if_pdb_path = pdb_file
            print('saving interface structures')
            create_pymol_session(pdb_filepath1=pathogen_if_pdb_path, shared_prot_chain1=row['pathogen_ppi_bait_chain'],
                                 interaction_id1=row['pathogen_interaction_id'], pdb_filepath2=human_if_pdb_path,
                                 shared_prot_chain2=row['human_ppi_bait_chain'], interaction_id2=row['human_interaction_id'],
                                 output_filepath=interface_output_filepath)
            os.remove(pathogen_if_pdb_path)
            os.remove(human_if_pdb_path)

    avg_plddt_df = pd.DataFrame(data=avg_plddt_dict)
    avg_plddt_df.to_csv(os.path.join(output_dir, 'avg_if_plddt.csv'), index=False)

    ss_df = pd.DataFrame(data=ss_dict)
    ss_df['pathogen_percentage_none'] = (ss_df['pathogen_-'] / ss_df['num_pathogen_if_residues']) * 100
    ss_df['pathogen_most_frequent_ss'] = ss_df[['pathogen_H', 'pathogen_B', 'pathogen_E', 'pathogen_G', 'pathogen_I', 'pathogen_T', 'pathogen_S', 'pathogen_-']].idxmax(axis=1)
    ss_df['human_percentage_none'] = (ss_df['human_-'] / ss_df['num_human_if_residues']) * 100
    ss_df['human_most_frequent_ss'] = ss_df[
        ['human_H', 'human_B', 'human_E', 'human_G', 'human_I', 'human_T', 'human_S', 'human_-']].idxmax(axis=1)

    ss_df.to_csv(os.path.join(output_dir, 'interfaces_secondary_struct.csv'), index=False)


if __name__ == '__main__':
    analyze_top_scores(output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis2/high_interface_comparison_scores',
                       only_top_human_interactor=True)