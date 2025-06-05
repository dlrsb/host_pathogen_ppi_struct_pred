import os
import inspect
import numpy as np
import pandas as pd


def interaction_id(uniprot_id_1, uniprot_id_2):
    # sorted pair of uniprot_id-s to compare interacting pairs from different sources
    return '_'.join(sorted([uniprot_id_1, uniprot_id_2]))


def workpath(path):
    # dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)


def uf(x):
    return '{:,}'.format(x)


def printlen(x, *args, **kwargs):
    name_ = inspect.stack()[1][3] #https://stackoverflow.com/questions/5067604/determine-function-name-from-within-that-function-without-using-traceback
    if name_ != '<module>':
        print(f'{name_}:', uf(len(x)), *args, **kwargs)
    else:
        print(uf(len(x)), *args, **kwargs)


def read_summary_source(output_filepath):
    df = pd.read_csv(workpath('23.10.02_dburke_kcl_OneDrive/summary_source.out.bz2'),
                      na_values={'pdockq_fd': 'chain'}, nrows=None).rename({'#protdir': 'protdir'}, axis=1)
    printlen(df, 'raw models')

    with open(workpath(f'23.10.02_dburke_kcl_OneDrive/af2-models-split.txt'), 'r') as f:
        af2_models_list = [line.strip() for line in f.readlines()]

    def parse_(s_): # this will discard non-dimers
        l_id_ = s_.split('/')[1].split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'

    def get_pdb_path(pair_id):
        pth1 = workpath(f'23.10.02_dburke_kcl_OneDrive/af2-models-split/{pair_id.split("_")[0]}/{pair_id}.pdb')
        pair_id2 = f'{pair_id.split("_")[1]}_{pair_id.split("_")[0]}'
        pth2 = workpath(f'23.10.02_dburke_kcl_OneDrive/af2-models-split/{pair_id.split("_")[1]}/{pair_id2}.pdb')
        if f'af2-models-split/{pair_id.split("_")[0]}/{pair_id}.pdb' in af2_models_list:
            return pth1
        elif f'af2-models-split/{pair_id.split("_")[1]}/{pair_id2}.pdb' in af2_models_list:
            return pth2
        else:
            return np.nan

    df.insert(0, 'interaction_id', [* df['protdir'].map(parse_) ])
    df.insert(2, 'folding_method', df['protdir'].map(lambda s: s.split('/')[-1]))
    df = df.query('interaction_id != "."').copy()
    df[['uniprot_id_A', 'uniprot_id_B']] = df['interaction_id'].str.split('_', expand=True)
    print(uf(len(df)), '\t', uf(len(set(df['interaction_id']))), 'models (interactions) after discarding non-dimers')
    df = df.drop_duplicates(subset=['interaction_id'], keep='first')
    print(df.shape)
    # df['pdb'] = df.apply(
    #     lambda r: workpath(f'23.10.02_dburke_kcl_OneDrive/af2-models-split/{r.interaction_id.split("_")[0]}/{r.interaction_id}.pdb'), axis=1)
    df['pdb'] = df['interaction_id'].apply(get_pdb_path)
    print(df.shape)
    df.dropna(axis=0, subset=['pdb'], inplace=True)
    print(df.shape)
    df.to_csv(output_filepath, index=False)
    return df


def read_ppi_reselect(output_filepath):
    fp_ = workpath('23.12.06_ppi_reselect/interface_best2_p10.csv') #8A
    # fp_ = workpath('23.12.06_ppi_reselect/interface_relaxed_best2_p10.csv') # 10A
    #fp_ = workpath('23.12.06_ppi_reselect/interface_strict_best2_p10.csv')  # 5A
    df = pd.read_csv(fp_, delim_whitespace=True).rename({
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

    df.insert(0, 'interaction_id', [*df['pair'].map(parse_)])
    df = df.query('interaction_id != "."').copy()
    df[['protein1', 'protein2']] = df['pair'].str.split('_', expand=True)
    df = df.drop_duplicates(subset=['interaction_id'], keep='first')
    df['pdb'] = df.apply(
        lambda r: workpath(f'23.12.06_ppi_reselect/af2-models-split/{r.pair.split("_")[0]}/{r.pair}.pdb'), axis=1)
    df.to_csv(output_filepath, index=False)
    return df.reset_index(drop=True)


def add_info(data_filepath, uniprot_info_filepath, output_filepath):
    data_df = pd.read_csv(data_filepath)
    uniprot_df = pd.read_table(uniprot_info_filepath, header=0, sep='\t')
    uniprot_df = uniprot_df[~uniprot_df['From'].duplicated(keep='first')]
    uniprot_df['families_domains'] = uniprot_df[
        ['CDD', 'Gene3D', 'InterPro', 'PROSITE', 'Pfam', 'NCBIfam', 'SMART', 'SUPFAM', 'PRINTS', 'HAMAP']].fillna('').agg(
        ';'.join, axis=1)
    organisms_dict = dict(uniprot_df[['From', 'Organism']].values)
    protein_names_dict = dict(uniprot_df[['From', 'Protein names']].values)
    domains_dict = dict(uniprot_df[['From', 'families_domains']].values)

    def get_species(uniprot_id):
        if uniprot_id in uniprot_df['From'].unique().tolist():
            search_id = uniprot_id
        else:
            search_id = uniprot_id.split('-')[0]
        try:
            species = organisms_dict[search_id].split(' (')[0]
            if 'subsp.' in species:
                species = species.split(' subsp.')[0]
        except Exception as e:
            print(e)
            print(f'Failed: {uniprot_id}\t{search_id}')
            species = 'NA'
        return species

    def get_names(uniprot_id):
        if uniprot_id in uniprot_df['From'].unique().tolist():
            search_id = uniprot_id
        else:
            search_id = uniprot_id.split('-')[0]
        try:
            names = protein_names_dict[search_id]
        except Exception as e:
            print(e)
            print(uniprot_id)
            names = 'NA'
        return names

    def get_domains(uniprot_id):
        if uniprot_id in uniprot_df['From'].unique().tolist():
            search_id = uniprot_id
        else:
            search_id = uniprot_id.split('-')[0]
        try:
            domains = domains_dict[search_id]
        except Exception as e:
            print(e)
            print(uniprot_id)
            domains = 'NA'
        return domains

    data_df[['organism1', 'organism2']] = data_df[['protein1', 'protein2']].applymap(get_species)
    data_df[['protein_name1', 'protein_name2']] = data_df[['protein1', 'protein2']].applymap(get_names)
    data_df[['families_domains1', 'families_domains2']] = data_df[['protein1', 'protein2']].applymap(get_domains)
    data_df.to_csv(output_filepath, index=False)


def filter_and_label_interactions(ppi_results_filepath, output_filepath):
    df = pd.read_csv(ppi_results_filepath)
    df_selected = df[(df['organism1'] == 'Homo sapiens') & (df['organism2'] == 'Homo sapiens')]
    df_selected['interaction_type'] = 'human-human'
    df_selected.to_csv(output_filepath, index=False)


if __name__ == '__main__':
    #read_summary_source(output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/other_organisms_af2_models.csv')
    read_ppi_reselect(output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/other_organisms_af2_models_ppi_reselect.csv')
    add_info(data_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect.csv',
             uniprot_info_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/idmapping_organisms_names_domains.tsv',
             output_filepath='/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info.csv')
    filter_and_label_interactions('/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info.csv',
                                  '/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/af2_other_species/other_organisms_af2_models_ppi_reselect_uniprot_info_human-human.csv')