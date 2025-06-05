import os
import copy

import pandas as pd
import tqdm, tqdm.contrib.concurrent

from shared_interface_analysis import SharedInterfaceAnalysis, parse_resid
from utils import calc_avg_plddt, get_interface_secondary_struct, create_pymol_session

OUTPUT_DIR_PATH = ''


def slurm_ntasks():
    return int(os.environ['SLURM_NTASKS'])


def prep_krogan_lab_results(results_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    print(results_df.shape)

    def get_model_filepath(x):
        model_path = 'unrelaxed_{}.pdb'.format(x['ranked0_model'])
        pth = '/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/{}/{}/{}/{}'.format(x['dataset'], x['name'], x['name'], model_path)
        return pth

    results_df['ranked0_model_path'] = results_df.apply(get_model_filepath, axis=1)

    results_df['pathogen_protein'] = results_df['protein1']
    results_df['human_protein'] = results_df['protein2'] # in these datasets, protein2 is always the human protein

    # Add chains - chains in unrelaxed PDB files always start with 'B' instead of 'A'
    results_df['chain1'] = 'B'
    results_df['chain2'] = 'C'
    results_df['interaction_type'] = 'host-pathogen'
    print(results_df.shape)

    results_df.to_csv(output_filepath, index=False)


def prep_hpidb_results(results_filepath, metadata_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    print(results_df.shape)
    metadata_df = pd.read_csv(metadata_filepath)
    # df_metadata['id'] = df_metadata['id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace('NP_', 'NP-')
    # df_metadata['pathogen_protein_id'] = df_metadata['pathogen_protein_id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace('NP_', 'NP-')
    # df_metadata['human_protein_id'] = df_metadata['human_protein_id'].str.replace('INTACT:', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace('NP_', 'NP-')
    # df_results['name'] = df_results['name'].str.replace('INTACT_', 'INTACT-').str.replace('PRO_', 'PRO-').str.replace('NP_', 'NP-')
    results_df['interaction_type'] = 'host-pathogen'
    results_df['chain1'] = 'A'
    results_df['chain2'] = 'B'
    results_df = results_df.merge(metadata_df[['id', 'human_protein_id', 'pathogen_protein_id']], how='inner', left_on='name', right_on='id')
    results_df['protein1'] = results_df['human_protein_id']
    results_df['protein2'] = results_df['pathogen_protein_id']
    results_df.rename(columns={'human_protein_id': 'human_protein', 'pathogen_protein_id': 'pathogen_protein'}, inplace=True)
    print(results_df.shape)
    results_df.to_csv(output_filepath, index=False)


def concat_host_pathogen_results(results_filepath1, results_filepath2, output_filepath):
    df1 = pd.read_csv(results_filepath1)
    df2 = pd.read_csv(results_filepath2)
    print(df1.shape)
    print(df2.shape)
    concat_df = pd.concat([df1, df2], axis=0)
    print(concat_df.shape)
    concat_df.to_csv(output_filepath, index=False)


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
            #print('Running Foldseek ComplexSearch for interfaces (all vs all)')
            #shared_interface_analysis_obj.calc_foldseek_interface_all_vs_all(workdir=foldseek_workdir)
            #shutil.rmtree(foldseek_workdir)
            if shared_interface_analysis_obj.df_interactors.shape[0] > 2:
                print('Clustering interfaces')
                shared_interface_analysis_obj.linkage(t=.9)
                shared_interface_analysis_obj.clustermap(os.path.join(OUTPUT_DIR_PATH, f'interface_clusters_{shared_interface_analysis_obj.uniprot_id_bait}.svg'))
                print('Creating PyMol session')
                shared_interface_analysis_obj.to_pymol(os.path.join(OUTPUT_DIR_PATH, f'interface_clusters_{shared_interface_analysis_obj.uniprot_id_bait}.pse'))
            results_df = shared_interface_analysis_obj.df_interface_scores
        else:
            results_df = pd.DataFrame()
    except Exception as e:
        print(f'Failed bait {shared_interface_analysis_obj.uniprot_id_bait}')
        print(e)
        results_df = pd.DataFrame()

    return results_df


def analyze_interfaces(host_pathogen_ppi_results_filepath, human_human_models_filepath, output_dir,
                       output_filename='apms_interface_comparison_scores.csv', iptm_threshold=0.5):
    global OUTPUT_DIR_PATH
    OUTPUT_DIR_PATH = output_dir

    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
    hp_ppi_data_df_selected = copy.deepcopy(hp_ppi_data_df[hp_ppi_data_df['iptm'] >= iptm_threshold])
    hp_ppi_data_df_selected.rename({'name': 'interaction_id',
                                    'huintaf2_pdockq_8Acutoff': 'pdockq',
                                    'huintaf2_interface_residues_chain1_8Acutoff': 'residues1',
                                    'huintaf2_interface_residues_chain2_8Acutoff': 'residues2',
                                    'ranked0_model_path': 'pdb'}, axis=1, inplace=True)
    hp_ppi_data_df_selected['interactor_protein_origin'] = 'pathogen'

    human_human_df = pd.read_csv(human_human_models_filepath)
    human_human_df_selected = human_human_df.query('pdockq >= 0.49') # CAPRI class "Medium" or "High" quality models only

    host_protein_ids = hp_ppi_data_df_selected['human_protein'].unique().tolist()

    df_models = pd.concat([hp_ppi_data_df_selected, human_human_df_selected], axis=0)

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

    bait_objs = [SharedInterfaceAnalysis(df_interactors_selected[df_interactors_selected['bait_id'] == uniprot_id]) for uniprot_id in host_protein_ids]

    ialign_dir='/cluster/scratch/dbaptista/ialign_tmpdir'
    if not os.path.exists(ialign_dir):
        os.mkdir(ialign_dir)

    all_scores = tqdm.contrib.concurrent.process_map(calc_scores, bait_objs, max_workers=slurm_ntasks(), chunksize=1)
    all_scores_df = pd.concat(all_scores, axis=0, ignore_index=True)
    all_scores_df.to_csv(os.path.join(output_dir, output_filename), index=False)


def label_interactors(interactor_id, interactors_df):
    interaction_type = interactors_df.loc[interactors_df['interactor_id'] == interactor_id, 'interaction_type'].unique().tolist()[0]
    if interaction_type == 'host-pathogen':
        return 'pathogen'
    else:
        return 'human'


def get_pathogen_human_proteins(row):
    if row['interactor1_origin'] == 'pathogen':
        pathogen_prot = row['interactor1']
        human_prot = row['interactor2']
    else:
        pathogen_prot = row['interactor2']
        human_prot = row['interactor1']
    return pathogen_prot, human_prot


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

    avg_plddt_df = pd.DataFrame(data=avg_plddt_dict)
    interactors_df_plddt = interactors_df.merge(avg_plddt_df, how='left', on=['bait_id', 'interactor_id'])

    ss_df = pd.DataFrame(data=ss_dict)
    ss_df['ss_percentage_none'] = (ss_df['-'] / ss_df['num_if_residues']) * 100
    ss_df['most_frequent_ss'] = ss_df[['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']].idxmax(axis=1)
    interactors_df_plddt_ss = interactors_df_plddt.merge(ss_df, how='left', on=['bait_id', 'interactor_id'])
    interactors_df_plddt_ss.to_csv(output_filepath, index=False)

    return interactors_df_plddt_ss


def analyze_pathogen_human_interface_scores(interface_scores_filepath, interactors_filepath, output_dir,
                                            only_top_human_interactor=False, metric_select_top='ialign8_itmscore'):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plddt_dirs = ['very_low_if_plddt', 'low_if_plddt', 'high_if_plddt', 'very_high_if_plddt']
    for plddt_dir in plddt_dirs:
        full_path = os.path.join(output_dir, plddt_dir)
        if not os.path.exists(full_path):
            os.mkdir(full_path)

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

        # Save PyMol session with full structures
        if not os.path.exists(full_struct_output_filepath):
            print('saving full structures')
            print(row['pathogen_protein'])
            print(row['human_protein'])
            print(row['pathogen_ppi_pdb'])
            print(row['pathogen_ppi_bait_chain']),
            print(row['pathogen_ppi_interaction_id']),
            print(row['human_ppi_pdb'])
            print(row['human_ppi_bait_chain'])
            print(row['human_ppi_interaction_id'])
            print(full_struct_output_filepath)
            create_pymol_session(pdb_filepath1=row['pathogen_ppi_pdb'], shared_prot_chain1=row['pathogen_ppi_bait_chain'],
                                 interaction_id1=row['pathogen_ppi_interaction_id'], pdb_filepath2=row['human_ppi_pdb'],
                                 shared_prot_chain2=row['human_ppi_bait_chain'], interaction_id2=row['human_ppi_interaction_id'],
                                 output_filepath=full_struct_output_filepath)


def pathogens_interface_similarity(pathogen_pathogen_interface_scores_filepath, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    scores_df = pd.read_csv(pathogen_pathogen_interface_scores_filepath)
    jaccard_above_030 = copy.deepcopy(scores_df[scores_df['jaccard_index'] >= 0.3])
    jaccard_above_030['pair1'] = jaccard_above_030['interactor1'] + '_' + jaccard_above_030['bait']
    jaccard_above_030['pair2'] = jaccard_above_030['interactor2'] + '_' + jaccard_above_030['bait']
    jaccard_above_030.to_csv(os.path.join(output_dir, 'pathogen_interactor_pairs_jaccard_030.csv'), index=False)
    pairs_with_similar_interfaces = list(set(jaccard_above_030['pair1'].to_list() + jaccard_above_030['pair2'].to_list() ))
    with open(os.path.join(output_dir, 'pathogen_interface_similarity.txt'), 'w') as f:
        f.write(str(len(pairs_with_similar_interfaces))+ '\n')


def calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath, output_dir, virus_only=False, iptm_cutoff=None,
                                    t=0.9):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    hp_ppi_data_df = pd.read_csv(host_pathogen_ppi_results_filepath)
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
            bait_ = SharedInterfaceAnalysis(df_interactors_selected.query('bait_id == @protein_id'))
            if bait_.df_interactors.shape[0] > 1:
                bait_.build_matrix()
                print('Calculating Jaccard Index')
                bait_.calc_jaccard_similarity()
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


def num_pathogen_proteins_per_bait(interactors_filepath, output_filepath):
    df = pd.read_csv(interactors_filepath)

    counts = {'bait':[], 'count': []}
    baits = df['bait_id'].unique().tolist()
    for bait in baits:
        df_selected = df[(df['bait_id'] == bait) & (df['interaction_type'] == 'host-pathogen')]
        counts['bait'].append(bait)
        counts['count'].append(df_selected.shape[0])

    counts_df = pd.DataFrame(data=counts)
    counts_df.sort_values(by='count', ascending=False, axis=0, inplace=True, ignore_index=True)
    counts_df.to_csv(output_filepath, index=False)


def get_novel_pathogen_interfaces(results_filepath, interactors_filepath, output_filepath, maximum_jaccard=0):
    df = pd.read_csv(results_filepath)
    interactors_df = pd.read_csv(interactors_filepath)
    baits = df['bait'].unique().tolist()
    novel_if_dict = {'bait': [], 'interactor': []}
    for bait in baits:
        df_selected = df[df['bait'] == bait]
        interactors = list(set(df_selected['interactor1'].to_list() + df_selected['interactor2'].to_list()))
        for interactor in interactors:
            df_interactor = df_selected[(df_selected['interactor1'] == interactor) | (df_selected['interactor2'] == interactor)]
            if df_interactor['jaccard_index'].le(maximum_jaccard).all():
                print(bait)
                print(df_interactor['jaccard_index'])
                print(df_interactor['jaccard_index'].unique().shape)
                novel_if_dict['bait'].append(bait)
                novel_if_dict['interactor'].append(interactor)

    novel_if_df = pd.DataFrame(data=novel_if_dict)
    novel_if_df = novel_if_df.merge(interactors_df[['bait_id', 'interactor_id', 'interaction_type']], how='left', left_on=['bait', 'interactor'], right_on=['bait_id', 'interactor_id'])
    novel_if_df.drop(columns=['bait_id', 'interactor_id'], inplace=True)
    novel_pathogen_if = novel_if_df[novel_if_df['interaction_type'] == 'host-pathogen']
    novel_pathogen_if.to_csv(output_filepath, index=False)


def count_interface_similarity_type(classified_interfaces_filepath):
    df = pd.read_csv(classified_interfaces_filepath)

    df_same = df[(df['similarity_type'] == 'same domain family') & (df['pathogen_interactor_avg_if_plddt'] > 50)]
    print(df_same.shape)
    print(f'Same domain family: {df_same["sorted_hp_pair"].unique().shape[0]} interactions')


    pairs_same = df_same['sorted_hp_pair'].unique().tolist()

    df_diff = df[(df['similarity_type'] == 'different domain, similar structural motif') & (
            df['pathogen_interactor_avg_if_plddt'] >= 50)]
    print(df_diff.shape)
    df_diff_sel = df_diff[~df_diff['sorted_hp_pair'].isin(pairs_same)]
    print(df_diff_sel.shape)
    print(f'Different domain, similar structural motif: {df_diff_sel["sorted_hp_pair"].unique().shape[0]} interactions')

    pairs_diff = df_diff_sel['sorted_hp_pair'].unique().tolist()
    pairs_same_diff = pairs_same + pairs_diff

    df_linear = df[df['similarity_type'].str.contains('linear')]
    print(df_linear.shape)
    df_linear_sel = df_linear[~df_linear['sorted_hp_pair'].isin(pairs_same_diff)]
    print(df_linear_sel.shape)
    print(f'Linear/less structured: {df_linear_sel["sorted_hp_pair"].unique().shape[0]} interactions')


if __name__ == '__main__':
    prep_krogan_lab_results(results_filepath='../results/krogan_lab_datasets/only_best_model/all/results_published_datasets_onlybestmodel_deduplicated.csv',
                            output_filepath='../results/krogan_lab_datasets/only_best_model/all/results_published_datasets_onlybestmodel_deduplicated_with_info.csv')

    prep_hpidb_results(results_filepath='../results/hpidb/hpidb_filtered_af2_results_deduplicated.csv',
                       metadata_filepath='../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered_deduplicated.csv',
                       output_filepath='../results/hpidb/hpidb_filtered_af2_results_deduplicated_with_info.csv')

    concat_host_pathogen_results(results_filepath1='../results/krogan_lab_datasets/only_best_model/all/results_published_datasets_onlybestmodel_deduplicated_with_info.csv',
                                 results_filepath2='../results/hpidb/hpidb_filtered_af2_results_deduplicated_with_info.csv',
                                 output_filepath='../results/merged_results/krogan_lab_hpidb_results_deduplicated.csv')

    analyze_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/merged_results/krogan_lab_hpidb_results_deduplicated.csv',
                       human_human_models_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_human-human.csv',
                       output_dir='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs',
                       output_filename='apms_interface_comparison_scores_allscores.csv',
                       iptm_threshold=0.5)

    add_interface_info(interactors_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected.csv',
                       output_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv')

    select_pathogen_vs_human_interactors(interface_scores_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_allscores.csv',
                                         interactors_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                                         output_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_human.csv')

    select_pathogen_vs_pathogen_interactors(interface_scores_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_allscores.csv',
                                            interactors_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                                            output_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_pathogen.csv')

    pathogens_interface_similarity('/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_pathogen.csv',
                                   output_dir='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/pathogen_pathogen')

    calc_unique_pathogen_interfaces(host_pathogen_ppi_results_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/merged_results/krogan_lab_hpidb_results_deduplicated.csv',
                                    output_dir='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/unique_pathogen_interfaces',
                                    iptm_cutoff=0.5,
                                    t=0.4)

    analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_human.csv',
                                            interactors_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                                            output_dir='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/top_jaccard',
                                            only_top_human_interactor=True,
                                            metric_select_top='jaccard_index')

    analyze_pathogen_human_interface_scores(interface_scores_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_human.csv',
                                            interactors_filepath='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                                            output_dir='/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/top_ialign',
                                            only_top_human_interactor=True,
                                            metric_select_top='ialign8_itmscore')

    num_pathogen_proteins_per_bait('/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                                   '/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/num_pathogen_proteins_per_bait.csv')

    get_novel_pathogen_interfaces('/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_allscores.csv',
                         '/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                         '/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/novel_pathogen_interfaces_maxjaccard0.csv')

    get_novel_pathogen_interfaces('/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/apms_interface_comparison_scores_pathogen_human.csv',
                         '/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/df_interactors_selected_with_interface_info.csv',
                         '/cluster/scratch/dbaptista/interface_analysis_deduplicated_pairs/pathogen_interfaces_not_targeted_by_human_prots_maxjaccard0.csv')

    count_interface_similarity_type('../results/interface_analysis_deduplicated_pairs/pathogen_human/interface_comparison_top_ialign_jaccard_classified.csv')