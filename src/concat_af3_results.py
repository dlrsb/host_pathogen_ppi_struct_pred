import pandas as pd


def prep_krogan_lab_results(results_filepath, output_filepath):
    results_df = pd.read_csv(results_filepath)
    print(results_df.columns)
    results_df['protein1'] = results_df['name'].str.split('_').str[0]
    results_df['protein2'] = results_df['name'].str.split('_').str[1]

    results_df['pathogen_protein'] = results_df['protein1']
    results_df['human_protein'] = results_df['protein2'] # in these datasets, protein2 is always the human protein

    # Add chains
    results_df['chain1'] = 'A'
    results_df['chain2'] = 'B'
    results_df['interaction_type'] = 'host-pathogen'
    results_df['dataset'] = 'AP-MS'

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
    results_df['dataset'] = 'HPIDB'
    results_df.to_csv(output_filepath, index=False)


def concat_host_pathogen_results(results_filepath1, results_filepath2, output_filepath):
    df1 = pd.read_csv(results_filepath1)
    df2 = pd.read_csv(results_filepath2)
    print(df1.shape)
    print(df2.shape)
    concat_df = pd.concat([df1, df2], axis=0)
    print(concat_df.shape)
    concat_df.to_csv(output_filepath, index=False)


if __name__ == '__main__':
    prep_krogan_lab_results(results_filepath="../results/krogan_lab_datasets/af3/krogan_lab_datasets_af3_results_deduplicated.csv",
                            output_filepath="../results/krogan_lab_datasets/af3/krogan_lab_datasets_af3_results_deduplicated_with_info.csv")
    prep_hpidb_results(results_filepath="../results/hpidb/af3/hpidb_af3_results_deduplicated.csv",
                       metadata_filepath="../data/hpidb2/hpidb_filtered_20250324/hpidb_selected_pairs_seqs_filtered_deduplicated.csv",
                       output_filepath="../results/hpidb/af3/hpidb_af3_results_deduplicated_with_info.csv")
    concat_host_pathogen_results(results_filepath1="../results/krogan_lab_datasets/af3/krogan_lab_datasets_af3_results_deduplicated_with_info.csv",
                                 results_filepath2="../results/hpidb/af3/hpidb_af3_results_deduplicated_with_info.csv",
                                 output_filepath="../results/merged_results/af3/krogan_lab_hpidb_results_deduplicated_af3.csv")
