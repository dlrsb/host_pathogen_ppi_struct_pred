import os
import glob
import copy
import numpy as np
from scipy.spatial.distance import pdist
import pandas as pd
from Bio import SeqIO

from shared_interface_analysis import SharedInterfaceAnalysis
from utils import load_structure
from plot_scores import classify_antibody


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


def get_human_ppi_models(uniprot_id, human_models_path):
    pdb_files = glob.glob(os.path.join(human_models_path, f'*{uniprot_id}*'))
    pdb_filepaths = [os.path.join(human_models_path, pdb_file) for pdb_file in pdb_files]
    return pdb_filepaths


def interaction_id(uniprot_id_1, uniprot_id_2):
    # sorted pair of uniprot_id-s to compare interacting pairs from different sources
    return '_'.join(sorted([uniprot_id_1, uniprot_id_2]))


def workpath(path):
    # dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)


def read_structures():
    #return pd.read_csv("C:\\Users\\dlrba\\Desktop\\structures_23.11.1.tsv", sep='\t')
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


def calc_percentage_shared_interface_old(interactor_dataframes, output_dir):
    host_pathogen_pair_counts = 0
    hp_with_shared_interface_count = 0
    for df in interactor_dataframes:
        host_pathogen_pair_counts += df['interaction_type'].value_counts()['host-pathogen']
        host_pathogen_pairs = df.loc[df['interaction_type'] == 'host-pathogen', 'interaction_id'].to_list()
        for hp_pair in host_pathogen_pairs:
            # for each pair with interaction_type == 'host-pathogen', check if its cluster ('label')
            # also has pairs with interaction_type = 'human-human'
            label = df.loc[df['interaction_id'] == hp_pair, 'labels'].item()
            if df[df['labels'] == label].value_counts().shape[0] > 1:
                hp_with_shared_interface_count += 1

    perc_same = (hp_with_shared_interface_count / host_pathogen_pair_counts) * 100
    perc_different = 100 - perc_same

    with open(os.path.join(output_dir, 'shared_interfaces_stats.txt'), 'w') as f:
        f.write(f'Total number of host-pathogen pairs considered:\t{host_pathogen_pair_counts}\n')
        f.write(f'Number of host-pathogen pairs which use interfaces similar to those used by human pairs:\t{hp_with_shared_interface_count}\n')
        f.write(f'Percentage of host-pathogen pairs which use interfaces similar to those used by human pairs:\t{perc_same}\n')
        f.write(f'Percentage of host-pathogen pairs which use different interfaces:\t{perc_different}\n')


def calc_percentage_shared_interface(jaccard_dataframes, pathogens_list, jaccard_threshold, output_dir):
    num_pathogen_pairs = 0
    num_pathogen_pairs_above_threshold = 0
    for df in jaccard_dataframes:
        pathogen_proteins = set(df['interactor1'].to_list() + df['interactor2'].to_list()).intersection(set(pathogens_list))
        num_pathogen_pairs += len(pathogen_proteins)
        for prot in pathogen_proteins:
            df_prot = df[df.eq(prot).any(axis=1)]
            print(df_prot)
            print(df_prot[df_prot['jaccard_index'] >= jaccard_threshold])
            if df_prot[df_prot['jaccard_index'] >= jaccard_threshold].shape[0] != 0:
                num_pathogen_pairs_above_threshold += 1

    perc_same = np.round((num_pathogen_pairs_above_threshold / num_pathogen_pairs) * 100, 2)

    with open(os.path.join(output_dir, f'shared_interfaces_stats_threshold_{jaccard_threshold}.txt'), 'w') as f:
        f.write(f'Total number of host-pathogen pairs considered:\t{num_pathogen_pairs}\n')
        f.write(f'Number of host-pathogen pairs which use interfaces similar to those used by human pairs (at least one human interactor with a Jaccard index over {jaccard_threshold}:\t{num_pathogen_pairs_above_threshold}\n')
        f.write(f'Percentage of host-pathogen pairs which use interfaces similar to those used by human pairs:\t{perc_same}%\n')


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
    #df_interactors.to_csv(os.path.join(output_dir, 'df_interactors.csv'), index=False)
    df_interactors_selected = df_interactors[df_interactors['bait_id'].isin(host_protein_ids)]
    df_interactors_selected.to_csv(os.path.join(output_dir, 'df_interactors_selected.csv'), index=False)

    jaccard_dfs = []
    for uniprot_id in host_protein_ids:
        try:
            bait_ = SharedInterfaceAnalysis(df_interactors_selected.query('bait_id == @uniprot_id'))
            #bait_.df_interactors.to_csv(os.path.join(output_dir, f'{uniprot_id}_df_interactors.csv'), index=False)
            #if bait_.df_matrix.shape[0] > 1:
            if bait_.df_interactors[bait_.df_interactors['interaction_type'] == 'human-human'].shape[0] >= 2:
                bait_.build_matrix()
                print('calculating jaccard index')
                jaccard_similarity = bait_.calc_jaccard_similarity(os.path.join(output_dir, f'{uniprot_id}_jaccard_similarity.csv'))
                jaccard_dfs.append(jaccard_similarity)
                bait_.calc_jaccard_similarity(os.path.join(output_dir, f'{uniprot_id}_jaccard_similarity_pathogens.csv'),
                                             only_compare_pathogen_against_human=True)
                print('calculating interface comparison scores')
                bait_.calc_interface_comparison_scores()
                bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_interface_comparison_scores.csv'))
                #bait_.calc_tmscore_interface(radius=5, workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/')
                #bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_interface_tmscores.csv'))
                #print('calculating TM-scores between interactors')
                #bait_.calc_tmscore(fname=os.path.join(output_dir, f'{uniprot_id}_pathogen_struct_alignment.csv'),
                #                   workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/')
                print('Running Foldseek ComplexSearch')
                bait_.calc_foldseek_interface(workdir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/foldseek_runs')
                bait_.save_scores_to_file(os.path.join(output_dir, f'{uniprot_id}_foldseek_scores.csv'))
                print('clustering interfaces')
                bait_.linkage(t=.9)
                bait_.clustermap(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.svg'))
                print('creating PyMol session')
                bait_.to_pymol(os.path.join(output_dir, f'interface_clusters_{uniprot_id}.pse'))
        except Exception as e:
            print(f'Failed bait {uniprot_id}')
            print(e)

    #Quantify % of cases where pathogens use the same interfaces as human proteins
    #print('calculating % of pathogens that use interfaces that are similar to the ones used by human interactors')
    #calc_percentage_shared_interface(jaccard_dfs, pathogen_protein_ids, 0.7, output_dir)
    #calc_percentage_shared_interface(jaccard_dfs, pathogen_protein_ids, 0.6, output_dir)
    #calc_percentage_shared_interface(jaccard_dfs, pathogen_protein_ids, 0.5, output_dir)


if __name__ == '__main__':
    # prep_host_pathogen_results(results_dir='../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model',
    #                            output_file='../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    # analyze_interfaces('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
    #                    output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/shared_interface_analysis2',
    #                    iptm_threshold=0.5)
    #prep_benchmarks_results('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/benchmarks_nonantibody_onlybestmodel.csv')
    analyze_interfaces('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv',
                       benchmark_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/benchmarks_nonantibody_onlybestmodel.csv',
                       output_dir='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/interface_analysis2',
                       iptm_threshold=0.5)
