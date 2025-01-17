import os
import glob
import copy
import shutil

import numpy as np
import pandas as pd
from Bio import SeqIO
import tqdm, tqdm.contrib.concurrent

from shared_interface_analysis import SharedInterfaceAnalysis, SharedInterfaceAnalysis2, parse_resid
from utils import load_structure, calc_avg_plddt, get_interface_secondary_struct, get_ialign_interface, create_pymol_session
from plot_scores import classify_antibody

OUTPUT_DIR_PATH = ''


def slurm_ntasks():
    return int(os.environ['SLURM_NTASKS'])


def prep_host_pathogen_results(results_dir, output_file):
    results_files = glob.glob(os.path.join(results_dir, '*.csv'))
    if any('results_published_datasets_onlybestmodel.csv' in x for x in results_files):
        results_files.remove(os.path.join(results_dir, 'results_published_datasets_onlybestmodel.csv'))
    if any('results_all_onlybestmodel.csv' in x for x in results_files):
        results_files.remove(os.path.join(results_dir, 'results_all_onlybestmodel.csv'))

    def get_filepath(x):
        model_path = 'unrelaxed_{}.pdb'.format(x['ranked0_model'])
        pth = '/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/{}/{}/{}/{}'.format(x['dataset'], x['name'], x['name'], model_path)
        return pth

    dataframes = []
    exclude = ['EV71-HEK293T', 'HBV-HUH7', 'IAV-HTBE']
    for res_file in results_files:
        if 'sars_mers_sars2_human_ppi' not in res_file:
            dataset = os.path.split(res_file)[-1].split('_')[0]
        else:
            dataset = 'sars_mers_sars2_human_ppi'

        if dataset not in exclude:
            df = pd.read_csv(res_file)
            df['dataset'] = dataset
            df['ranked0_model_path'] = df.apply(get_filepath, axis=1)
            dataframes.append(df)
    concat_df = pd.concat(dataframes, axis=0)
    # Remove results for models with nan coordinates (empty AF2 scores columns):
    concat_df = concat_df.dropna(subset=['ptm', 'iptm', 'ranking_confidence'], axis=0)
    # Drop duplicate protein pairs (protein pairs that were tested in more than 1 cell line):
    #concat_df = concat_df.drop_duplicates(subset='name', keep='first')
    concat_df = concat_df.sort_values('iptm', ascending=False).drop_duplicates(subset='name', keep='first')
    # Add chains - chains in unrelaxed PDB files always start with 'B' instead of 'A'
    concat_df['chain1'] = 'B'
    concat_df['chain2'] = 'C'
    concat_df.to_csv(output_file, index=False)


def prep_benchmarks_results(output_file):

    def get_filepath(x):
        model_path = 'unrelaxed_{}_chains_renamed.pdb'.format(x['ranked0_model'])
        pth = '/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/{}/{}/{}/{}'.format(x['dataset'], x['name'], x['name'], model_path)
        return pth

    def get_chains(x):
        struct = load_structure(x['ranked0_model_path'])
        chains_list = [chain.id for chain in struct.get_chains()]
        return chains_list

    def get_uniprot_ids(x):
        data_dir = '_'.join(x['dataset'].split('_')[0:-1])
        fasta_pth = f'/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/{data_dir}/fullseq_fasta_files/protein_pairs/{x["assembly_id"]}.fasta'
        seqrecords = list(SeqIO.parse(fasta_pth, 'fasta'))
        uniprot1 = ''
        uniprot2 = ''
        for seqrecord in seqrecords:
            description_split = seqrecord.description.split(' ')
            if description_split[0] == x['chain1']:
                uniprot1 = description_split[-1]
            elif description_split[0] == x['chain2']:
                uniprot2 = description_split[-1]
        return [uniprot1, uniprot2]

    def find_pathogen_protein(x):
        data_dir = '_'.join(x['dataset'].split('_')[0:-1])
        fasta_pth = f'/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/{data_dir}/fullseq_fasta_files/protein_pairs/{x["assembly_id"]}.fasta'
        seqrecords = list(SeqIO.parse(fasta_pth, 'fasta'))
        pathogen_protein = ''
        human_protein = ''
        for sr in seqrecords:
            description_split = sr.description.split(' ')
            uniprot_cols = x[x == description_split[-1]].index.to_list()
            print(uniprot_cols)
            if len(uniprot_cols) == 1: # this happens when the FASTA description does not have a Uniprot ID
                pathogen_protein = description_split[-1]
            else:
                for col in uniprot_cols:
                    if 'uniprots.uniprot_id' in col:
                        entity_id = col.split('.')[1]
                        organism_col = f'polymer_entity_instances.{entity_id}.polymer_entity.uniprots.uniprot_source_organism_scientific_name'
                        if x[organism_col] == 'Homo sapiens': # if organism is not Homo sapiens, this chain is the pathogen protein
                            human_protein = description_split[-1]
                        else:
                            pathogen_protein = description_split[-1]

        return human_protein, pathogen_protein

    bacteria_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv')
    bacteria_df['dataset'] = 'bacteria_mammalia_dimers_fullseq'
    bacteria_metadata_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/70identity_search_results.csv',)
    bacteria_df_with_metadata = bacteria_df.merge(bacteria_metadata_df, how='inner', left_on='name', right_on='assembly_id')
    virus_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv')
    virus_df['dataset'] = 'virus_mammalia_dimers_fullseq'
    virus_metadata_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/70identity_search_results.csv',)
    virus_df_with_metadata = virus_df.merge(virus_metadata_df, how='inner', left_on='name', right_on='assembly_id')
    concat_df = pd.concat([bacteria_df_with_metadata, virus_df_with_metadata], axis=0, ignore_index=True)
    concat_df['ranked0_model_path'] = concat_df.apply(get_filepath, axis=1)
    concat_df[['chain1', 'chain2']] = concat_df.apply(get_chains, axis=1, result_type='expand')
    concat_df[['protein1', 'protein2']] = concat_df.apply(get_uniprot_ids, axis=1, result_type='expand')
    concat_df['antibody'] = concat_df.apply(classify_antibody, axis=1)
    concat_df['interaction_type'] = 'host-pathogen'
    # Select only entries where the host protein is from Homo sapiens:
    selected_df = concat_df[(concat_df['antibody'] == 'Not antibody') & ((concat_df['polymer_entity_instances.1.polymer_entity.uniprots.uniprot_source_organism_scientific_name'] == 'Homo sapiens') | (concat_df['polymer_entity_instances.2.polymer_entity.uniprots.uniprot_source_organism_scientific_name'] == 'Homo sapiens'))]
    selected_df[['human_protein', 'pathogen_protein']] = selected_df.apply(find_pathogen_protein, axis=1, result_type='expand')
    selected_df.to_csv(output_file, index=False)


def interaction_id(uniprot_id_1, uniprot_id_2):
    # sorted pair of uniprot_id-s to compare interacting pairs from different sources
    return '_'.join(sorted([uniprot_id_1, uniprot_id_2]))


def workpath(path):
    # dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)


def read_structures():
    return pd.read_csv(workpath('23.11.01_human_protein_map/structures_23.11.1.tsv'), sep='\t')


def af2_uniprot_id():
    # AF2 human single-fragment structures (set used for pocket detection, interface modelling, etc)
    return set(read_structures()['uniprot_id'])


def read_ppi_reselect():
    fp_ = workpath('23.12.06_ppi_reselect/interface_best2_p10.csv') #8A
    # fp_ = workpath('23.12.06_ppi_reselect/interface_relaxed_best2_p10.csv') # 10A
    #fp_ = workpath('23.12.06_ppi_reselect/interface_strict_best2_p10.csv')  # 5A
    df_ = pd.read_csv(fp_, delim_whitespace=True).rename({
        'chain': 'chain1',
        'residues': 'residues1',
        'chain': 'chain1',
        'chain.1': 'chain2',
        'residues.1': 'residues2',
    }, axis=1)

    def parse_(s_):
        l_id_ = s_.split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'

    df_.insert(0, 'interaction_id', [*df_['pair'].map(parse_)])
    df_ = df_.query('interaction_id != "."').copy()
    df_[['pair1', 'pair2']] = df_['pair'].str.split('_', expand=True)
    # return df_
    df_ = df_.query('(protein1 in @af2_uniprot_id()) & (protein2 in @af2_uniprot_id())').copy()
    df_ = df_.drop_duplicates(subset=['interaction_id'], keep='first')
    df_['pdb'] = df_.apply(
        lambda r: workpath(f'23.12.06_ppi_reselect/af2-models-split/{r.pair.split("_")[0]}/{r.pair}.pdb'), axis=1)
    return df_.reset_index(drop=True)


def analyze_interfaces(host_pathogen_ppi_results_filepath, benchmark_results_filepath, output_dir, iptm_threshold=0.5,
                       use_benchmark_datasets=False):
    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
    hp_ppi_data_df['interaction_type'] = 'host-pathogen'
    hp_ppi_data_df_selected = copy.deepcopy(hp_ppi_data_df[hp_ppi_data_df['iptm'] >= iptm_threshold])
    hp_ppi_data_df_selected[['protein1', 'protein2']] = hp_ppi_data_df_selected['name'].str.split(pat='_', expand=True)
    hp_ppi_data_df_selected.rename({'name': 'interaction_id',
                                    'huintaf2_pdockq_8Acutoff': 'pdockq',
                                    'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                    'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                    'ranked0_model_path': 'pdb'}, axis=1, inplace=True)

    human_ppi_data_df = read_ppi_reselect().query('pdockq > .5')
    human_ppi_data_df['interaction_type'] = 'human-human'
    #human_ppi_data_df.to_csv(os.path.join(output_dir, 'human_ppi_data_df.csv'), index=False)

    if use_benchmark_datasets:
        benchmark_data_df = pd.read_csv(benchmark_results_filepath)
        benchmark_data_df_selected = copy.deepcopy(benchmark_data_df[benchmark_data_df['dockq'] >= 0.23]) #or should I use a different metric/threshold??
        benchmark_data_df_selected.rename({'name': 'interaction_id',
                                           'huintaf2_pdockq_8Acutoff': 'pdockq',
                                           'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                           'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                           'ranked0_model_path': 'pdb'}, axis=1, inplace=True)

        host_protein_ids = list(set(hp_ppi_data_df_selected['protein2'].to_list() + benchmark_data_df_selected['human_protein'].to_list()))
        pathogen_protein_ids = list(set(hp_ppi_data_df_selected['protein1'].to_list() + benchmark_data_df_selected['pathogen_protein'].to_list()))
        df_models = pd.concat([hp_ppi_data_df_selected, benchmark_data_df_selected, human_ppi_data_df], axis=0, ignore_index=True)
    else:
        host_protein_ids = hp_ppi_data_df_selected['protein2'].unique().tolist() # protein2 is always the host protein
        pathogen_protein_ids = hp_ppi_data_df_selected['protein1'].unique().tolist()
        df_models = pd.concat([hp_ppi_data_df_selected, human_ppi_data_df], axis=0)

    cols_ = ['bait_id', 'bait_ifresid', 'bait_chain', 'interactor_id', 'interactor_chain', 'interaction_id', 'pdockq',
             'pdb', 'interaction_type', 'interactor_ifresid']
    q_ne_ = 'protein1 != protein2'
    q_eq_ = 'protein1 == protein2'
    df_interactors = pd.concat([
        df_models.query(q_ne_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_eq_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_ne_).rename(
            {'protein2': 'bait_id', 'residues2': 'bait_ifresid', 'chain2': 'bait_chain', 'chain1': 'interactor_chain',
             'protein1': 'interactor_id', 'residues1': 'interactor_ifresid'},
            axis=1)[cols_],
    ], axis=0)
    df_interactors_selected = df_interactors[df_interactors['bait_id'].isin(host_protein_ids)]
    df_interactors_selected.to_csv(os.path.join(output_dir, 'df_interactors_selected_new.csv'), index=False)

    for uniprot_id in host_protein_ids:
        if not os.path.exists(os.path.join(output_dir, f'{uniprot_id}_jaccard_similarity.csv')):
            try:
                bait_ = SharedInterfaceAnalysis(df_interactors_selected.query('bait_id == @uniprot_id'))
                if bait_.df_interactors.shape[0] > 1:
                    bait_.build_matrix()
                    #print('calculating jaccard index')
                    bait_.calc_jaccard_similarity(os.path.join(output_dir, f'{uniprot_id}_jaccard_similarity.csv'))
                    #bait_.calc_jaccard_similarity(os.path.join(output_dir, f'{uniprot_id}_jaccard_similarity_pathogens.csv'),
                    #                             only_compare_pathogen_against_human=True)
                    #print('calculating interface comparison scores')
                    #bait_.calc_interface_comparison_scores()
                    #bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_interface_comparison_scores.csv'))
                    #bait_.calc_tmscore_interface(radius=5, workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/')
                    #bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_interface_tmscores.csv'))
                    #print('calculating TM-scores between interactors')
                    #bait_.calc_tmscore(fname=os.path.join(output_dir, f'{uniprot_id}_pathogen_struct_alignment.csv'),
                    #                   workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/')

                    #print('Running Foldseek ComplexSearch')
                    #bait_.calc_foldseek_interface(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/foldseek_runs')
                    #bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_foldseek_scores.csv'))

                    if bait_.df_interactors[bait_.df_interactors['interaction_type'] == 'human-human'].shape[0] >= 2:
                        print('clustering interfaces')
                        bait_.linkage(t=.9)
                        bait_.clustermap(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.svg'))
                        print('creating PyMol session')
                        bait_.to_pymol(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.pse'))
            except Exception as e:
                print(f'Failed bait {uniprot_id}')
                print(e)


def analyze_interfaces2(host_pathogen_ppi_results_filepath, output_dir,
                        output_filename='apms_interface_comparison_scores_testfoldseek.csv', iptm_threshold=0.5):
    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
    hp_ppi_data_df['interaction_type'] = 'host-pathogen'
    hp_ppi_data_df_selected = copy.deepcopy(hp_ppi_data_df[hp_ppi_data_df['iptm'] >= iptm_threshold])
    hp_ppi_data_df_selected[['protein1', 'protein2']] = hp_ppi_data_df_selected['name'].str.split(pat='_', expand=True)
    hp_ppi_data_df_selected.rename({'name': 'interaction_id',
                                    'huintaf2_pdockq_8Acutoff': 'pdockq',
                                    'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                    'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                    'ranked0_model_path': 'pdb'}, axis=1, inplace=True)
    hp_ppi_data_df_selected['interactor_protein_origin'] = 'pathogen'

    human_ppi_data_df = read_ppi_reselect().query('pdockq > .5') # TODO: change this to iptm once we get these values from David
    human_ppi_data_df['interaction_type'] = 'human-human'

    #host_protein_ids = hp_ppi_data_df_selected['protein2'].unique().tolist() # protein2 is always the host protein
    host_protein_ids = ['P24941']
    df_models = pd.concat([hp_ppi_data_df_selected, human_ppi_data_df], axis=0)

    cols_ = ['bait_id', 'bait_ifresid', 'bait_chain', 'interactor_id', 'interactor_chain', 'interaction_id', 'pdockq',
             'pdb', 'interactor_ifresid', 'interaction_type']
    q_ne_ = 'protein1 != protein2'
    q_eq_ = 'protein1 == protein2'
    df_interactors = pd.concat([
        df_models.query(q_ne_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_eq_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_ne_).rename(
            {'protein2': 'bait_id', 'residues2': 'bait_ifresid', 'chain2': 'bait_chain', 'chain1': 'interactor_chain',
             'protein1': 'interactor_id', 'residues1': 'interactor_ifresid'},
            axis=1)[cols_],
    ], axis=0)

    df_interactors_selected = df_interactors[df_interactors['bait_id'].isin(host_protein_ids)]
    df_interactors_selected.to_csv(os.path.join(output_dir, 'df_interactors_selected.csv'), index=False)

    all_scores = []

    for uniprot_id in host_protein_ids:
        try:
            bait_ = SharedInterfaceAnalysis2(df_interactors_selected.query('bait_id == @uniprot_id'))
            if bait_.df_interactors.shape[0] > 1:
                #bait_.build_matrix()
                #print('Calculating jaccard index')
                #bait_.calc_jaccard_similarity()
                #print('Calculating iAlign scores')
                #bait_.calc_ialign_scores(scoring_metric='tm', distance_cutoff=8)
                #bait_.calc_ialign_scores_all_vs_all(scoring_metric='tm', distance_cutoff=8)
                # print('Calculating interface TM-score')
                # bait_.calc_tmscore(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/',
                #                    interface_only=True, radius=5)
                print('Running Foldseek ComplexSearch for interfaces (pairwise)')
                bait_.calc_foldseek_interface(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/foldseek_runs')
                # print('Running Foldseek ComplexSearch for interfaces (all vs all)')
                # bait_.calc_foldseek_interface_all_vs_all(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/foldseek_runs',
                #                                          foldseek_output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/foldseek_output.csv')
                all_scores.append(bait_.df_interface_scores)

                # if bait_.df_interactors[bait_.df_interactors['interaction_type'] == 'human-human'].shape[0] >= 2:
                #     print('Clustering interfaces')
                #     bait_.linkage(t=.9)
                #     bait_.clustermap(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.svg'))
                #     print('Creating PyMol session')
                #     bait_.to_pymol(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.pse'))

        except Exception as e:
            print(f'Failed bait {uniprot_id}')
            print(e)

    all_scores_df = pd.concat(all_scores, axis=0, ignore_index=True)
    all_scores_df.to_csv(os.path.join(output_dir, output_filename), index=False)


def calc_scores(shared_interface_analysis_obj):
    try:
        if shared_interface_analysis_obj.df_interactors.shape[0] > 1:
            shared_interface_analysis_obj.build_matrix()
            print('Calculating jaccard index')
            shared_interface_analysis_obj.calc_jaccard_similarity()
            print('Calculating iAlign scores')
            shared_interface_analysis_obj.calc_ialign_scores_all_vs_all(scoring_metric='tm', distance_cutoff=8,
                                                                        ialign_dir='/cluster/scratch/dbaptista/ialign_tmpdir')
            shared_interface_analysis_obj.calc_ialign_scores_all_vs_all(scoring_metric='is', distance_cutoff=8,
                                                                        ialign_dir='/cluster/scratch/dbaptista/ialign_tmpdir')
            # print('Calculating interface TM-score')
            # shared_interface_analysis_obj.calc_tmscore(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/',
            #                    interface_only=True, radius=5)
            #foldseek_workdir = f'/cluster/scratch/dbaptista/foldseek_pathogen_if/{shared_interface_analysis_obj.uniprot_id_bait}'
            #os.mkdir(foldseek_workdir)
            # print('Running Foldseek ComplexSearch for interfaces (pairwise)')
            # shared_interface_analysis_obj.calc_foldseek_interface(workdir=foldseek_workdir)
            #print('Running Foldseek ComplexSearch for interfaces (all vs all)')
            #shared_interface_analysis_obj.calc_foldseek_interface_all_vs_all(workdir=foldseek_workdir)
            #shutil.rmtree(foldseek_workdir)
            # if shared_interface_analysis_obj.df_interactors[shared_interface_analysis_obj.df_interactors['interaction_type'] == 'human-human'].shape[0] >= 2:
            #     print('Clustering interfaces')
            #     shared_interface_analysis_obj.linkage(t=.9)
            #     shared_interface_analysis_obj.clustermap(os.path.join(OUTPUT_DIR_PATH, f'interface_clusters_{shared_interface_analysis_obj.uniprot_id_bait}.svg'))
            #     print('Creating PyMol session')
            #     shared_interface_analysis_obj.to_pymol(os.path.join(OUTPUT_DIR_PATH, f'interface_clusters_{shared_interface_analysis_obj.uniprot_id_bait}.pse'))
            results_df = shared_interface_analysis_obj.df_interface_scores
        else:
            results_df = pd.DataFrame()
    except Exception as e:
        print(f'Failed bait {shared_interface_analysis_obj.uniprot_id_bait}')
        print(e)
        results_df = pd.DataFrame()

    return results_df


def filter_and_label_interactions(ppi_results_filepath, output_filepath):
    df = pd.read_csv(ppi_results_filepath)
    df_selected = df[(df['organism1'] == 'Homo sapiens') & (df['organism2'] == 'Homo sapiens')]
    df_selected['interaction_type'] = 'human-human'

    # OLD:
    # def label_row(row):
    #     if row['organism1'] == 'Homo sapiens' and row['organism2'] == 'Homo sapiens':
    #         interaction_type = 'human-human'
    #     elif row['organism1'] == 'Homo sapiens' and row['organism2'] != 'Homo sapiens':
    #         interaction_type = 'host-pathogen'
    #     elif row['organism1'] != 'Homo sapiens' and row['organism2'] == 'Homo sapiens':
    #         interaction_type = 'host-pathogen'
    #     else:
    #         interaction_type = 'other'
    #     return interaction_type

    # df['interaction_type'] = df.apply(label_row, axis=1)
    # df2 = df[~(df['interaction_type'] == 'other')] # excludes pairs that don't have at least 1 human protein
    # exclude = ['Citrobacter rodentium', 'Woodchuck hepatitis virus', 'Oryctolagus cuniculus',
    #            'Schizosaccharomyces pombe', 'Murid herpesvirus 4', 'Aotus azarae', 'Canis lupus familiaris',
    #            'Trichosanthes kirilowii', 'Simian immunodeficiency virus', 'Arabidopsis thaliana',
    #            'Caenorhabditis elegans', 'Mus musculus', 'Rhipicephalus sanguineus', 'Chlamydia caviae',
    #            'Chlorocebus sabaeus', 'Rattus norvegicus', 'Drosophila melanogaster', 'Chicken anemia virus',
    #            'Bovine papillomavirus type 1', 'Sus scrofa', 'Gallus gallus', 'Chlorocebus aethiops', 'Ovis aries',
    #            'Lens culinaris', 'Equine herpesvirus 2', 'Xenopus laevis', 'Cricetulus griseus',
    #            'Rhesus monkey rhadinovirus H26-95', 'Saccharomyces cerevisiae', 'Bos taurus', np.nan]
    # exclude pairs with proteins from organisms that are not human pathogens:
    #df_selected = df2[~((df2['organism1'].isin(exclude)) | (df2['organism2'].isin(exclude)))]
    df_selected.to_csv(output_filepath, index=False)


def analyze_interfaces3(host_pathogen_ppi_results_filepath, other_models_results_filepath, output_dir,
                        output_filename='apms_interface_comparison_scores.csv', iptm_threshold=0.5):
    global OUTPUT_DIR_PATH
    OUTPUT_DIR_PATH = output_dir

    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
    hp_ppi_data_df['interaction_type'] = 'host-pathogen'
    hp_ppi_data_df_selected = copy.deepcopy(hp_ppi_data_df[hp_ppi_data_df['iptm'] >= iptm_threshold])
    hp_ppi_data_df_selected[['protein1', 'protein2']] = hp_ppi_data_df_selected['name'].str.split(pat='_', expand=True)
    hp_ppi_data_df_selected.rename({'name': 'interaction_id',
                                    'huintaf2_pdockq_8Acutoff': 'pdockq',
                                    'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                    'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                    'ranked0_model_path': 'pdb'}, axis=1, inplace=True)
    hp_ppi_data_df_selected['interactor_protein_origin'] = 'pathogen'

    #human_ppi_data_df = read_ppi_reselect().query('pdockq > .5') # TODO: change this to iptm once we get these values from David
    #human_ppi_data_df['interaction_type'] = 'human-human'
    other_models_df = pd.read_csv(other_models_results_filepath)
    other_models_df_selected = other_models_df.query('pdockq >= 0.49')


    host_protein_ids = hp_ppi_data_df_selected['protein2'].unique().tolist() # protein2 is always the host protein
    # other_models_hp_pairs_df = other_models_df_selected[other_models_df_selected['organism1'] != other_models_df_selected['organism2']]
    # for i, r in other_models_hp_pairs_df.iterrows():
    #     if r['organism1'] == 'Homo sapiens':
    #         host_protein_ids.append(r['protein1'])
    #     elif r['organism2'] == 'Homo sapiens':
    #         host_protein_ids.append(r['protein2'])
    # host_protein_ids = list(set(host_protein_ids))

    #df_models = pd.concat([hp_ppi_data_df_selected, human_ppi_data_df], axis=0)
    df_models = pd.concat([hp_ppi_data_df_selected, other_models_df_selected], axis=0)

    cols_ = ['bait_id', 'bait_ifresid', 'bait_chain', 'interactor_id', 'interactor_chain', 'interaction_id', 'pdockq',
             'pdb', 'interactor_ifresid', 'interaction_type']
    q_ne_ = 'protein1 != protein2'
    q_eq_ = 'protein1 == protein2'
    df_interactors = pd.concat([
        df_models.query(q_ne_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_eq_).rename(
            {'protein1': 'bait_id', 'residues1': 'bait_ifresid', 'chain1': 'bait_chain', 'chain2': 'interactor_chain',
             'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
            axis=1)[cols_],
        df_models.query(q_ne_).rename(
            {'protein2': 'bait_id', 'residues2': 'bait_ifresid', 'chain2': 'bait_chain', 'chain1': 'interactor_chain',
             'protein1': 'interactor_id', 'residues1': 'interactor_ifresid'},
            axis=1)[cols_],
    ], axis=0)

    df_interactors_selected = df_interactors[df_interactors['bait_id'].isin(host_protein_ids)]
    df_interactors_selected.to_csv(os.path.join(output_dir, 'df_interactors_selected.csv'), index=False)
    print(df_interactors_selected)

    bait_objs = [SharedInterfaceAnalysis2(df_interactors_selected[df_interactors_selected['bait_id'] == uniprot_id]) for uniprot_id in host_protein_ids]

    all_scores = tqdm.contrib.concurrent.process_map(calc_scores, bait_objs, max_workers=slurm_ntasks(), chunksize=1)
    all_scores_df = pd.concat(all_scores, axis=0, ignore_index=True)
    all_scores_df.to_csv(os.path.join(output_dir, output_filename), index=False)


def label_interactors(interactor_id, interactors_df):
    interaction_type = interactors_df.loc[interactors_df['interactor_id'] == interactor_id, 'interaction_type'].unique().tolist()[0]
    if interaction_type == 'host-pathogen':
        return 'pathogen'
    else:
        return 'human'


def select_pathogen_vs_human_interactors(interface_scores_filepath, interactors_filepath, output_filepath):
    scores_df = pd.read_csv(interface_scores_filepath)
    interactors_df = pd.read_csv(interactors_filepath)
    scores_df['interactor1_origin'] = scores_df['interactor1'].apply(label_interactors, args=(interactors_df,))
    scores_df['interactor2_origin'] = scores_df['interactor2'].apply(label_interactors, args=(interactors_df,))
    scores_df_pathogen_human = scores_df[scores_df['interactor1_origin'] != scores_df['interactor2_origin']]
    scores_df_pathogen_human.to_csv(output_filepath, index=False)
    return scores_df_pathogen_human


def select_pathogen_vs_pathogen_interactors(interface_scores_filepath, interactors_filepath, output_filepath):
    scores_df = pd.read_csv(interface_scores_filepath)
    interactors_df = pd.read_csv(interactors_filepath)
    scores_df['interactor1_origin'] = scores_df['interactor1'].apply(label_interactors, args=(interactors_df,))
    scores_df['interactor2_origin'] = scores_df['interactor2'].apply(label_interactors, args=(interactors_df,))
    scores_df_pathogen_pathogen = scores_df[(scores_df['interactor1_origin'] == 'pathogen') & (scores_df['interactor2_origin'] == 'pathogen')]
    scores_df_pathogen_pathogen.to_csv(output_filepath, index=False)
    return scores_df_pathogen_pathogen


def add_interface_info(interactors_filepath, output_filepath):
    interactors_df = pd.read_csv(interactors_filepath)

    avg_plddt_dict = {'bait_id': [], 'interactor_id': [], 'interactor_avg_if_plddt': []}
    ss_dict = {'bait_id': [], 'interactor_id': [], 'H': [], 'B': [], 'E': [], 'G': [], 'I': [], 'T': [], 'S': [], '-': [],
               'num_if_residues': []}

    for index, row in interactors_df.iterrows():
        avg_plddt = calc_avg_plddt(pdb_filepath=row['pdb'],
                                   interactor_chain=row['interactor_chain'],
                                   selected_residues=sorted(list(parse_resid(row['interactor_ifresid']))))
        print(avg_plddt)
        avg_plddt_dict['bait_id'].append(row['bait_id'])
        avg_plddt_dict['interactor_id'].append(row['interactor_id'])
        avg_plddt_dict['interactor_avg_if_plddt'].append(avg_plddt)

        ss_dict['bait_id'].append(row['bait_id'])
        ss_dict['interactor_id'].append(row['interactor_id'])
        ss_dict['num_if_residues'].append(len(parse_resid(row['interactor_ifresid'])))
        dssp = get_interface_secondary_struct(pdb_filepath=row['pdb'], chain=row['interactor_chain'],
                                              interface_residues=parse_resid(row['interactor_ifresid']))
        for key in dssp:
            ss_dict[key].append(dssp[key])

    print(avg_plddt_dict)
    avg_plddt_df = pd.DataFrame(data=avg_plddt_dict)
    interactors_df_plddt = interactors_df.merge(avg_plddt_df, how='left', on=['bait_id', 'interactor_id'])

    ss_df = pd.DataFrame(data=ss_dict)
    ss_df['ss_percentage_none'] = (ss_df['-'] / ss_df['num_if_residues']) * 100
    ss_df['most_frequent_ss'] = ss_df[['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']].idxmax(axis=1)
    interactors_df_plddt_ss = interactors_df_plddt.merge(ss_df, how='left', on=['bait_id', 'interactor_id'])
    interactors_df_plddt_ss.to_csv(output_filepath, index=False)
    print(interactors_df_plddt_ss.columns)

    return interactors_df_plddt_ss


def analyze_pathogen_human_interface_scores(interface_scores_filepath, interactors_filepath, output_dir,
                                            only_top_human_interactor=False, metric_select_top='ialign8_itmscore'):
    plddt_dirs = ['very_low_if_plddt', 'low_if_plddt', 'high_if_plddt', 'very_high_if_plddt']
    for plddt_dir in plddt_dirs:
        full_path = os.path.join(output_dir, plddt_dir)
        if not os.path.exists(full_path):
            os.mkdir(full_path)

    def get_pathogen_human_proteins(row):
        if row['interactor1_origin'] == 'pathogen':
            pathogen_prot = row['interactor1']
            human_prot = row['interactor2']
        else:
            pathogen_prot = row['interactor2']
            human_prot = row['interactor1']
        return pathogen_prot, human_prot

    scores_df = pd.read_csv(interface_scores_filepath)
    scores_df[['pathogen_protein', 'human_protein']] = scores_df.apply(get_pathogen_human_proteins, axis=1,
                                                                       result_type='expand')

    scores_df['hp_pair'] = scores_df['pathogen_protein'] + '_' + scores_df['bait']
    interactors_df = pd.read_csv(interactors_filepath)

    scores_df_selected = scores_df[(scores_df['ialign8_itmscore'] >= 0.5) & (scores_df['jaccard_index'] >= 0.3)]
    scores_df_selected.to_csv(os.path.join(output_dir, 'pathogen_human_interface_comparison_selected.csv'),
                              index=False)

    num_pairs_jaccard_030 = scores_df.loc[scores_df['jaccard_index'] >= 0.3, 'hp_pair'].unique().shape[0]
    num_pairs_selected = scores_df_selected['hp_pair'].unique().shape[0]

    with open(os.path.join(output_dir, 'pathogen_human_if_scores_summary.txt'), 'w') as f:
        f.write(f'# of host-pathogen protein pairs with Jaccard Index >= 0.3: {num_pairs_jaccard_030} \n')
        f.write(f'# of host-pathogen protein pairs with Jaccard Index >= 0.3 & iAlign iTM-score >= 0.5: {num_pairs_selected} \n')

    if only_top_human_interactor:
        df_grouped_max = scores_df_selected.loc[scores_df_selected.groupby(by=['pathogen_protein', 'bait'])[metric_select_top].idxmax()]
        scores_df_selected = df_grouped_max
        df_grouped_max.to_csv(os.path.join(output_dir, f'pathogen_human_interface_comparison_selected_only_top_{metric_select_top}.csv'),
                              index=False)

    merged_df1 = scores_df_selected.merge(interactors_df, how='left', left_on=['bait', 'pathogen_protein'],
                                   right_on=['bait_id', 'interactor_id'])
    merged_df1.rename({x: f'pathogen_ppi_{x}' for x in interactors_df.columns}, axis=1, inplace=True)
    merged_df2 = merged_df1.merge(interactors_df, how='left', left_on=['bait', 'human_protein'], right_on=['bait_id', 'interactor_id'])
    merged_df2.rename({x: f'human_ppi_{x}' for x in interactors_df.columns}, axis=1, inplace=True)

    print(merged_df2.columns)

    df_linear = merged_df2[(merged_df2['pathogen_ppi_most_frequent_ss'] == '-') & (merged_df2['human_ppi_most_frequent_ss'] == '-')]
    df_linear.to_csv(os.path.join(output_dir, 'pathogen_human_interface_comparison_linear_ss.csv'))

    for index, row in merged_df2.iterrows():
        print(row['pathogen_protein'])
        print(row['human_protein'])
        print(row['pathogen_ppi_pdb'])

        avg_plddt = row['pathogen_ppi_interactor_avg_if_plddt']

        if avg_plddt < 50.0:
            output_dir2 = os.path.join(output_dir, 'very_low_if_plddt')
        elif (avg_plddt >= 50.0) and (avg_plddt < 70.0):
            output_dir2 = os.path.join(output_dir, 'low_if_plddt')
        elif (avg_plddt >= 70.0) and (avg_plddt < 90.0):
            output_dir2 = os.path.join(output_dir, 'high_if_plddt')
        else:
            output_dir2 = os.path.join(output_dir, 'very_high_if_plddt')

        full_struct_output_filepath = os.path.join(output_dir2, f'full_{row["pathogen_protein"]}_{row["human_protein"]}_bait{row["bait"]}.pse')
        interface_output_filepath = os.path.join(output_dir2, f'interface_{row["pathogen_protein"]}_{row["human_protein"]}_bait{row["bait"]}.pse')

        # Save PyMol session with full structures
        if not os.path.exists(full_struct_output_filepath):
            print('saving full structures')
            create_pymol_session(pdb_filepath1=row['pathogen_ppi_pdb'], shared_prot_chain1=row['pathogen_ppi_bait_chain'],
                                 interaction_id1=row['pathogen_ppi_interaction_id'], pdb_filepath2=row['human_ppi_pdb'],
                                 shared_prot_chain2=row['human_ppi_bait_chain'], interaction_id2=row['human_ppi_interaction_id'],
                                 output_filepath=full_struct_output_filepath)

        # # Run ialign with selected metric and get interface PDB file (PDB files with _int, in "outputs" folder)
        # if not os.path.exists(interface_output_filepath):
        #     print('running ialign to get interface PDB files')
        #     interface_output_dir = '/cluster/scratch/dbaptista/ialign_interface_pdbs'
        #     get_ialign_interface(structure1_filepath=row['pathogen_ppi_pdb'], structure2_filepath=row['human_ppi_pdb'],
        #                          scoring_metric='tm', normalization_method='ave', distance_cutoff=8,
        #                          if_outputdir=interface_output_dir)
        #     interface_pdb_files = glob.glob(os.path.join(interface_output_dir, '*_int.pdb'))
        #     print(interface_pdb_files)
        #     for pdb_file in interface_pdb_files:
        #         if row['human_ppi_interaction_id'] in pdb_file:
        #             human_if_pdb_path = pdb_file
        #         else:
        #             pathogen_if_pdb_path = pdb_file
        #     print('saving interface structures')
        #     create_pymol_session(pdb_filepath1=pathogen_if_pdb_path, shared_prot_chain1=row['pathogen_ppi_bait_chain'],
        #                          interaction_id1=row['pathogen_ppi_interaction_id'], pdb_filepath2=human_if_pdb_path,
        #                          shared_prot_chain2=row['human_ppi_bait_chain'], interaction_id2=row['human_ppi_interaction_id'],
        #                          output_filepath=interface_output_filepath)
        #     os.remove(pathogen_if_pdb_path)
        #     os.remove(human_if_pdb_path)


def pathogens_interface_similarity(p126913264269877athogen_pathogen_interface_scores_filepath, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    scores_df = pd.read_csv(pathogen_pathogen_interface_scores_filepath)
    jaccard_above_030 = copy.deepcopy(scores_df[scores_df['jaccard_index'] >= 0.3])
    jaccard_above_030['pair1'] = jaccard_above_030['interactor1'] + '_' + jaccard_above_030['bait']
    jaccard_above_030['pair2'] = jaccard_above_030['interactor2'] + '_' + jaccard_above_030['bait']
    jaccard_above_030.to_csv(os.path.join(output_dir, 'pathogen_interactor_pairs_jaccard_030.csv'), index=False)
    pairs_with_similar_interfaces = list(set(jaccard_above_030['pair1'].to_list() + jaccard_above_030['pair2'].to_list() ))
    print(len(pairs_with_similar_interfaces))


def calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath, output_dir, virus_only=False, iptm_cutoff=None,
                                    t=0.9):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
    hp_ppi_data_df['interaction_type'] = 'host-pathogen'
    hp_ppi_data_df[['protein1', 'protein2']] = hp_ppi_data_df['name'].str.split(pat='_', expand=True)
    hp_ppi_data_df.rename({'name': 'interaction_id',
                                    'huintaf2_pdockq_8Acutoff': 'pdockq',
                                    'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                    'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                    'ranked0_model_path': 'pdb'}, axis=1, inplace=True)
    hp_ppi_data_df['interactor_protein_origin'] = 'pathogen'

    if iptm_cutoff is not None:
        hp_ppi_data_df = hp_ppi_data_df[hp_ppi_data_df['iptm'] >= iptm_cutoff]

    if virus_only:
        df_pathogen = copy.deepcopy(hp_ppi_data_df[~hp_ppi_data_df['interaction_id'].str.contains('Chlamydia|Mtb|MERS|SARS-CoV1|H1N1|H3N2', regex=True)])
    else:
        df_pathogen = copy.deepcopy(hp_ppi_data_df[~hp_ppi_data_df['interaction_id'].str.contains('MERS|SARS-CoV1|H1N1|H3N2', regex=True)])

    df_pathogen_nonempty_if = df_pathogen[df_pathogen['pdockq'] != '[]']

    pathogen_protein_ids = df_pathogen_nonempty_if['protein1'].unique().tolist() # use to determine unique pathogen interfaces

    cols = ['bait_id', 'bait_ifresid', 'bait_chain', 'interactor_id', 'interactor_chain', 'interaction_id', 'pdockq',
             'pdb', 'interactor_ifresid', 'interaction_type']
    df_interactors = df_pathogen_nonempty_if.rename({'protein1': 'bait_id', 'residues1': 'bait_ifresid',
                                                  'chain1': 'bait_chain', 'chain2': 'interactor_chain',
                                                  'protein2': 'interactor_id', 'residues2': 'interactor_ifresid'},
                                                 axis=1)[cols]
    df_interactors_selected = df_interactors[df_interactors['bait_id'].isin(pathogen_protein_ids)]

    results = {'pathogen_protein': [], f'unique_interface_count_t{t}': []}
    for protein_id in pathogen_protein_ids:
        results['pathogen_protein'].append(protein_id)
        try:
            bait_ = SharedInterfaceAnalysis2(df_interactors_selected.query('bait_id == @protein_id'))
            if bait_.df_interactors.shape[0] > 1:
                bait_.build_matrix()
                print('Calculating Jaccard Index')
                bait_.calc_jaccard_similarity() # TODO: I don't think this part is necessary. Remove
                print('Clustering interfaces')
                bait_.linkage(t=t)  # previously  0.9 # TODO: try to decide what the best value of t might be
                num_unique_interfaces = bait_.df_interactors['labels'].unique().shape[0]
                results[f'unique_interface_count_t{t}'].append(num_unique_interfaces)
            else:
                results[f'unique_interface_count_t{t}'].append(1)

        except Exception as e:
            print(f'Failed protein {protein_id}')
            print(e)

    results_df = pd.DataFrame(data=results)
    if virus_only:
        name = 'viral'
    else:
        name = 'pathogen'
    csv_filepath = f'num_unique_{name}_interfaces_per_{name}_prot.csv'
    results_df.to_csv(os.path.join(output_dir, csv_filepath), index=False)
    with open(os.path.join(output_dir, 'num_unique_interfaces.txt'), 'w') as f:
        f.write(f'Total number of {name}-human pairs considered: {df_pathogen["interaction_id"].unique().shape[0]}\n')
        f.write(f'Total number of {name}-human pairs where the {name} protein has a defined interface (8A, all-atom): {df_pathogen_nonempty_if["interaction_id"].unique().shape[0]}\n')
        f.write(f'Total number of {name} proteins: {df_pathogen_nonempty_if["protein1"].unique().shape[0]}\n')
        f.write(f'Number of unique {name} interfaces (t={t}, criterion=distance): {sum(results[f"unique_interface_count_t{t}"])}\n')


if __name__ == '__main__':
    # prep_host_pathogen_results(results_dir='../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model',
    #                            output_file='../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    #prep_benchmarks_results('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/benchmarks_nonantibody_onlybestmodel.csv')
    # analyze_interfaces('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                    benchmark_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/benchmarks_nonantibody_onlybestmodel.csv',
    #                    output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis2',
    #                    iptm_threshold=0.5,
    #                    use_benchmark_datasets=False)
    # analyze_interfaces2('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                     output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3',
    #                     output_filename='apms_interface_comparison_scores_testfoldseek_pairwise.csv',
    #                     iptm_threshold=0.5)
    # select_pathogen_vs_human_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/apms_interface_comparison_scores.csv',
    #                                      interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/df_interactors_selected.csv',
    #                                      output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/apms_interface_comparison_scores_pathogen_human.csv')
    # select_pathogen_vs_pathogen_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/apms_interface_comparison_scores.csv',
    #                                      interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/df_interactors_selected.csv',
    #                                      output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/apms_interface_comparison_scores_pathogen_pathogen.csv')
    # add_interface_info(interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/df_interactors_selected.csv',
    #                    output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/df_interactors_selected_with_interface_info2.csv')
    # # select_pathogen_vs_human_interactors(interface_scores_filepath="C:\\Users\\dlrba\\Desktop\\apms_interface_comparison_scores.csv",
    #                                      interactors_filepath="C:\\Users\\dlrba\\Desktop\\df_interactors_selected.csv",
    #                                      output_filepath="C:\\Users\\dlrba\\Desktop\\apms_interface_comparison_scores_pathogen_human.csv")
    # select_pathogen_vs_pathogen_interactors(interface_scores_filepath="C:\\Users\\dlrba\\Desktop\\apms_interface_comparison_scores.csv",
    #                                         interactors_filepath="C:\\Users\\dlrba\\Desktop\\df_interactors_selected.csv",
    #                                         output_filepath="C:\\Users\\dlrba\\Desktop\\apms_interface_comparison_scores_pathogen_pathogen.csv")
    # pathogens_interface_similarity("C:\\Users\\dlrba\\Desktop\\apms_interface_comparison_scores_pathogen_pathogen.csv",
    #                                output_dir="C:\\Users\\dlrba\\Desktop")
    # calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                                 output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/unique_viral_interfaces')
    # analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/apms_interface_comparison_scores_pathogen_human.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3/df_interactors_selected_with_interface_info.csv',
    #                                         output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis3',
    #                                         only_top_human_interactor=True)

    # analyze_interfaces3('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                     output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4',
    #                     output_filename='apms_interface_comparison_scores_allscores.csv',
    #                     iptm_threshold=0.5)
    # add_interface_info(interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/df_interactors_selected.csv',
    #                    output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/df_interactors_selected_with_interface_info.csv')
    # select_pathogen_vs_human_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_allscores.csv',
    #                                      interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/df_interactors_selected_with_interface_info.csv',
    #                                      output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_pathogen_human.csv')
    # select_pathogen_vs_pathogen_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_allscores.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/df_interactors_selected_with_interface_info.csv',
    #                                         output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_pathogen_pathogen.csv')
    # pathogens_interface_similarity('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_pathogen_human.csv',
    #                                output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/pathogen_pathogen')
    # calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                                 output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/unique_viral_interfaces')
    # analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/apms_interface_comparison_scores_pathogen_human.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4/df_interactors_selected_with_interface_info.csv',
    #                                         output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis4',
    #                                         only_top_human_interactor=True)

    # filter_and_label_interactions(ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info.csv',
    #                               output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_filtered.csv')
    # analyze_interfaces3(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                     other_models_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_filtered.csv',
    #                     output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5',
    #                     output_filename='apms_interface_comparison_scores_allscores.csv',
    #                     iptm_threshold=0.5)
    # add_interface_info(interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected.csv',
    #                    output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected_with_interface_info.csv')
    # select_pathogen_vs_human_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_allscores.csv',
    #                                      interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected_with_interface_info.csv',
    #                                      output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_pathogen_human.csv')
    # select_pathogen_vs_pathogen_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_allscores.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected_with_interface_info.csv',
    #                                         output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_pathogen_pathogen.csv')
    # pathogens_interface_similarity('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_pathogen_pathogen.csv',
    #                                output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/pathogen_pathogen')
    # calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                                 output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/unique_viral_interfaces') # TODO: might need to change to add in host-pathogen pairs in DB's dataset
    # analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_pathogen_human.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected_with_interface_info.csv',
    #                                         output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5',
    #                                         only_top_human_interactor=True)
    # calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                                 output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/unique_pathogen_interfaces',
    #                                 iptm_cutoff=0.5,
    #                                 t=0.4)

    # analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/apms_interface_comparison_scores_pathogen_human.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/df_interactors_selected_with_interface_info.csv',
    #                                         output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis5/top_jaccard',
    #                                         only_top_human_interactor=True,
    #                                         metric_select_top='jaccard_index')

    # NEW (all human-human models from David's dataset, no host-pathogen models:
    # filter_and_label_interactions(ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info.csv',
    #                               output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_human-human.csv')
    # analyze_interfaces3(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                     other_models_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_human-human.csv',
    #                     output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6',
    #                     output_filename='apms_interface_comparison_scores_allscores.csv',
    #                     iptm_threshold=0.5)
    # add_interface_info(interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected.csv',
    #                    output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected_with_interface_info.csv')
    # select_pathogen_vs_human_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_allscores.csv',
    #                                      interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected_with_interface_info.csv',
    #                                      output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_pathogen_human.csv')
    # select_pathogen_vs_pathogen_interactors(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_allscores.csv',
    #                                         interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected_with_interface_info.csv',
    #                                         output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_pathogen_pathogen.csv')
    # pathogens_interface_similarity('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_pathogen_pathogen.csv',
    #                                output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/pathogen_pathogen')
    # calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                                 output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/unique_pathogen_interfaces',
    #                                 iptm_cutoff=0.5,
    #                                 t=0.4)
    analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_pathogen_human.csv',
                                            interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected_with_interface_info.csv',
                                            output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/top_jaccard',
                                            only_top_human_interactor=True,
                                            metric_select_top='jaccard_index')
    analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/apms_interface_comparison_scores_pathogen_human.csv',
                                            interactors_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/df_interactors_selected_with_interface_info.csv',
                                            output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis6/top_ialign',
                                            only_top_human_interactor=True,
                                            metric_select_top='ialign8_itmscore')