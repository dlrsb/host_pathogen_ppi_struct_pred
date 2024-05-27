import ast
import copy
import itertools

import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn as sns


def calc_jaccard(residues1, residues2):
    return 1 - len(set(residues1).intersection(set(residues2))) / len(set(residues1).union(set(residues2)))


def calc_unique_viral_interfaces():
    #df = pd.read_csv('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/results_published_datasets_onlybestmodel.csv')
    df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    # df_iav = df[df['name'].str.contains('IAV')]
    # df_iav['subtype'] = df_iav['name'].str.split('_').str[0].str.split('-').str[-1]
    # df_iav['human_protein'] = df_iav['name'].str.split('_').str[-1]
    # h1n1_preys = set(df_iav.loc[df_iav['subtype'] == 'H1N1', 'human_protein'].to_list())
    # h5n1_preys = set(df_iav.loc[df_iav['subtype'] == 'H5N1', 'human_protein'].to_list())
    # h3n2_preys = set(df_iav.loc[df_iav['subtype'] == 'H3N2', 'human_protein'].to_list())
    # print(len(h1n1_preys))
    # print(len(h5n1_preys))
    # print(len(h3n2_preys))
    # print(len(h1n1_preys.union(h5n1_preys).union(h3n2_preys)))
    df_viral = copy.deepcopy(df[~df['name'].str.contains('Chlamydia|Mtb|MERS|SARS-CoV1|H1N1|H3N2', regex=True)])

    def get_pathogen_and_human_proteins(row):
        pathogen_protein, human_protein = row['name'].split('_')
        if 'IAV' in pathogen_protein:
            pathogen_protein_split = pathogen_protein.split('-')
            new_pathogen_protein_order = [pathogen_protein_split[0], pathogen_protein_split[-1]] + pathogen_protein_split[1:-1]
            pathogen_protein = '-'.join(new_pathogen_protein_order)
        # elif 'MERS' in pathogen_protein:
        #     pathogen_protein = '-'.join(['Betacoronavirus'] + pathogen_protein.split('-')[1:])
        # elif 'SARS' in pathogen_protein:
        #     pathogen_protein = '-'.join(['Betacoronavirus'] + pathogen_protein.split('-')[2:])
        return pathogen_protein, human_protein

    #df_viral[['pathogen_protein', 'human_protein']] = df_viral['name'].str.split('_', expand=True)
    df_viral[['pathogen_protein', 'human_protein']] = df_viral.apply(get_pathogen_and_human_proteins, axis=1, result_type='expand')
    df_viral_nonempty_if = df_viral[df_viral['huintaf2_interface_residues_chain1_8Acutoff'] != '[]']
    results = {'pathogen_protein': [], 'unique_interface_count': [], 'unique_interface_count_t0.9': [],
               'unique_interface_count_t0.4': []}
    for pathogen_prot in df_viral['pathogen_protein'].unique().tolist():
        df_pathogen = df_viral_nonempty_if[df_viral_nonempty_if['pathogen_protein'] == pathogen_prot]
        df_pathogen.to_csv(f"C:\\Users\\dlrba\\Desktop\\jaccard_viral_interfaces\\{pathogen_prot}_df.csv", index=False)
        results['pathogen_protein'].append(pathogen_prot)
        results['unique_interface_count'].append(df_pathogen["huintaf2_interface_residues_chain1_8Acutoff"].unique().shape[0])
        if df_pathogen.shape[0] == 1:
            results['unique_interface_count_t0.9'].append(1)
            results['unique_interface_count_t0.4'].append(1)
        else:
            pairwise_combinations = list(itertools.combinations(df_pathogen['name'].unique().tolist(), 2))
            jaccard_results = {'pair1': [], 'pair2': [], 'jaccard': []}
            for (pair1, pair2) in pairwise_combinations:
                jaccard_results['pair1'].append(pair1)
                jaccard_results['pair2'].append(pair2)
                res1 = ast.literal_eval(df_pathogen.loc[df_pathogen['name']==pair1, 'huintaf2_interface_residues_chain1_8Acutoff'].item())
                res2 = ast.literal_eval(df_pathogen.loc[df_pathogen['name']==pair2, 'huintaf2_interface_residues_chain1_8Acutoff'].item())
                jaccard = calc_jaccard(res1, res2)
                jaccard_results['jaccard'].append(jaccard)
            jaccard_results_df = pd.DataFrame(data=jaccard_results)
            jaccard_results_df.to_csv(f"C:\\Users\\dlrba\\Desktop\\jaccard_viral_interfaces\\{pathogen_prot}_jaccard.csv", index=False)
            print(pathogen_prot)
            #print(jaccard_results_df)
            linkage_ = scipy.cluster.hierarchy.linkage(jaccard_results_df['jaccard'].values, method='average')
            clusters_09 = scipy.cluster.hierarchy.fcluster(Z=linkage_, criterion='distance', t=0.9)
            clusters_04 = scipy.cluster.hierarchy.fcluster(Z=linkage_, criterion='distance', t=0.4)
            results['unique_interface_count_t0.9'].append(len(set(clusters_09)))
            results['unique_interface_count_t0.4'].append(len(set(clusters_04)))
            plt.figure()
            scipy.cluster.hierarchy.dendrogram(linkage_, color_threshold=0.9)
            plt.savefig(f"C:\\Users\\dlrba\\Desktop\\jaccard_viral_interfaces\\{pathogen_prot}_dendrogram.png")
            #plt.show()
            plt.close()
    results_df = pd.DataFrame(data=results)
    #print(results_df)
    print(f'Total number of viral-human pairs considered: {df_viral["name"].unique().shape[0]}')
    print(f'Total number of viral-human pairs where the viral protein has a defined interface (8A, all-atom): {df_viral_nonempty_if["name"].unique().shape[0]}')
    print(f'Total number of virus proteins: {df_viral_nonempty_if["pathogen_protein"].unique().shape[0]}')
    print(f'Number of unique viral interfaces: {sum(results["unique_interface_count"])}')
    print(f'Number of unique viral interfaces (t=0.9, criterion=distance): {sum(results["unique_interface_count_t0.9"])}')
    print(f'Number of unique viral interfaces (t=0.4, criterion=distance): {sum(results["unique_interface_count_t0.4"])}')
    results_df.to_csv('C:\\Users\\dlrba\\Desktop\\jaccard_viral_interfaces\\num_unique_viral_interfaces_per_viral_prot.csv', index=False)


def plot_num_unique_interfaces_per_virus(output_filepath):
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=14)  # legend fontsize
    plt.rc('figure', titlesize=18)  # fontsize of the figure title

    df = pd.read_csv("C:\\Users\\dlrba\\Desktop\\jaccard_viral_interfaces\\num_unique_viral_interfaces_per_viral_prot.csv")

    def get_pathogen(x):
        if 'IAV' in x:
            pathogen = '-'.join(x.split('-')[0:2])
        elif 'SARS-CoV2' in x:
            pathogen = 'SARS-CoV-2'
        elif 'ZIKVfp' in x:
            pathogen = 'ZIKV'
        else:
            pathogen = x.split('-')[0]
        return pathogen

    df['pathogen'] = df['pathogen_protein'].apply(get_pathogen)
    print(df['pathogen'])
    grouped = df.groupby(by='pathogen').sum().sort_values(by='unique_interface_count_t0.4', ascending=False)
    print(grouped.index)
    plt.figure(figsize=(5, 6))
    sns.barplot(x=grouped['unique_interface_count_t0.4'].to_list(), y=grouped.index.tolist(), orient='h', color='silver')
    plt.ylabel('Virus')
    plt.xlabel('Interfaces')
    #plt.show()
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    calc_unique_viral_interfaces()
    plot_num_unique_interfaces_per_virus(output_filepath='C:\\Users\\dlrba\\Desktop\\virus_num_interfaces.png')