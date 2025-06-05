import copy
import os
import datetime
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, mannwhitneyu, wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns

sns.color_palette('deep')
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.figure(figsize=(4, 4))


def get_capri_class(x):
    if 0.00 <= x < 0.23:
        capri_class = 'Incorrect'
    elif 0.23 <= x < 0.49:
        capri_class = 'Acceptable'
    elif 0.49 <= x < 0.80:
        capri_class = 'Medium'
    elif x >= 0.80:
        capri_class = 'High'
    return capri_class


def classify_release_date(release_date):
    date_str = release_date.split('T')[0]
    y, m, d = [int(x) for x in date_str.split('-')]
    date = datetime.date(y, m, d)
    cutoff_date = datetime.date(2021, 9, 30)
    if date > cutoff_date:
        return 'After training set cutoff'
    else:
        return 'Before training set cutoff'


def classify_antibody(row):
    description_cols = ['polymer_entity_instances.1.polymer_entity.rcsb_polymer_entity.pdbx_description',
                        'polymer_entity_instances.2.polymer_entity.rcsb_polymer_entity.pdbx_description']
    antibody_nanobody_terms = ['1-2C7', '3LRH intrabody', '5D', '7F', 'A20.1 VHH', 'A26.8 VHH', 'Ab08',
                               'Anti-F4+ETEC bacteria VHH variable region',
                               'Anti-Marburgvirus Nucleoprotein Single Domain Antibody A',
                               'Anti-Marburgvirus Nucleoprotein Single Domain Antibody C',
                               'Anti-Rev Antibody Fab single-chain variable fragment, light chain,Anti-Rev Antibody Fab single-chain variable fragment, heavy chain',
                               'Anti-Zaire ebolavirus Nucleoprotein Single Domain Antibody Zaire C (ZC)',
                               'Anti-Zaire ebolavirus Nucleoprotein Single Domain Antibody Zaire E (ZE)',
                               'Antibody C8 VHH domain', 'B39 VHH', 'C144 scFv', 'C5 nanobody',
                               'CA1698 camel antibody fragment', 'CAMELID VHH 9',
                               'Camelid heavy-chain antibody variable fragment cAb-F11N',
                               'Camelid heavy-chain antibody variable fragment cAb-G10S',
                               'Camelid heavy-chain antibody variable fragment cAb-H7S',
                               'DH270.6 single chain variable fragment',
                               'E3', 'F2 nanobody', 'H11-H4', 'Hi113 protein', 'Ig gamma-3 chain C region', 'JLE-E5',
                               'JLI-G10', 'JLI-H11', 'JLK-G12', 'JSG-C1', 'LaM8', 'Llama antibody D7', 'M2e-VHH-23m',
                               'Megabody 177', 'Monobody', 'Monobody S12', 'MraYAA nanobody', 'NB1A2', 'NB8194',
                               'NB_1B11', 'NanoE6', 'NanoF7', 'Nanobody', 'Nanobody (VHH) Nano-27',
                               'Nanobody (VHH) Nano-4', 'Nanobody 2', 'Nanobody 3-2A2-4', 'Nanobody 327', 'Nanobody 33',
                               'Nanobody 6', 'Nanobody 7', 'Nanobody 74', 'Nanobody 8', 'Nanobody B6', 'Nanobody DL4',
                               'Nanobody Nb-007', 'Nanobody Nb-ER14', 'Nanobody Nb-ER19', 'Nanobody Nb14527, NbSOS3',
                               'Nanobody NbFedF7', 'Nanobody P17', 'Nanobody VHH AA6', 'Nanobody VHH AH3',
                               'Nanobody VHH6', 'Nanobody9047', 'Nb113 Camel antibody fragment', 'Nb70', 'Nb8193',
                               'Nb_MsbA 1', 'Protein ca1697 (nanobody)', 'SCFV513',
                               'Single chain variable fragment of the non-neutralizing antibody DAO5',
                               'Single-chain Fv', 'Single-domain antibody 20ipaD', 'Single-domain antibody C2',
                               'VHH 1-21', 'VHH R303', 'VHH antiboby', 'VHH-12', 'VHH-C2', 'VHH6 nanobody',
                               'Variable domain of antibody scFv9', 'ZV-48 Antibody scFv',
                               'anti-HIV llama VHH antibody A12', 'anti-SARS-CoV-2 receptor binding domain VHH',
                               'anti-dengue Mab 4E11', 'anti-sars scFv antibody, 80R', 'antibody scFv',
                               'immunoglobulin heavy chain variable region', 'lama VHH antibody 2E7',
                               "mCherry's nanobody LaM4", 'n3113', 'nanobody', 'nanobody 282', 'nanobody NB4',
                               'nanobody SARS VHH-72', 'nb2b4', 'nb_1A7', 'neutralizing nanobody NM1226',
                               'neutralizing nanobody NM1230', 'scFv 2D10', 'scFv E4', 'scFv H2526',
                               'scFv of 9C12 antibody', 'single-chain Fv antibody fragment (scFv)',
                               'single-domain antibody JMK-E3', 'single-domain antibody JMK-H2']
    if any(row[description_cols[0]] == term for term in antibody_nanobody_terms) or any(row[description_cols[1]] == term for term in antibody_nanobody_terms):
        return 'Antibody'
    else:
        return 'Not antibody'


def plot_score_distributions(results_dfs_list, metric, output_filepath, hue_col=None, plot_type='hist_with_kde', bins=10,
                             figsize=(4,4)):
    # plot types available: 'hist', 'kde', 'hist_with_kde

    concat_df = pd.concat(results_dfs_list, axis=0, ignore_index=True).dropna(subset=metric)
    if plot_type == 'hist':
        g = sns.displot(data=concat_df, x=metric, bins=bins, hue=hue_col, height=4, aspect=1) # height 7.2
    elif plot_type == 'kde':
        g = sns.displot(data=concat_df, x=metric, kind='kde', hue=hue_col, height=4, aspect=1, clip=[0, 1])
    elif plot_type == 'hist_with_kde':
        g = sns.displot(data=concat_df, x=metric, kde=True, hue=hue_col, bins=bins, height=4, aspect=1)

    if metric == 'mmalign_tmscore':
        plt.xlabel('TM-score')
    elif metric == 'dockq':
        plt.xlabel('DockQ')

    ax = g.axes
    ax[0,0].set_yticks(np.arange(0, ax[0,0].get_ylim()[1], 25))

    if output_filepath is not None:
        plt.savefig(output_filepath, dpi=300,
                    bbox_inches='tight')
    plt.show()
    plt.close()


def plot_true_vs_predicted_scores(results_df, metric1, metric2, output_filepath,
                                  correlation_type='pearson', figsize=(4, 4), hue=None, color=None):
    metrics_df = results_df[[metric1, metric2]].dropna(axis=0)
    metric1_data = metrics_df[metric1]
    metric2_data = metrics_df[metric2]
    if correlation_type == 'pearson':
        correlation = pearsonr(metric1_data, metric2_data)[0]
    elif correlation_type == 'spearman':
        correlation = spearmanr(metric1_data, metric2_data)[0]

    # create scatterplot
    fig, ax = plt.subplots(figsize=figsize)
    if hue is not None:
        g = sns.scatterplot(data=results_df, x=metric1,
                            y=metric2, ax=ax, alpha=0.5, hue=hue)
        g.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    else:
        g = sns.scatterplot(x=metric1_data, y=metric2_data, ax=ax, alpha=0.5, color=color)

    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ',
                       'jaccard_index': 'Jaccard Index',
                       'ialign8_itmscore': 'iAlign iTM-score (8 Å cutoff)',
                       'iptm': 'ipTM',
                       'ptm': 'pTM',
                       'huintaf2_pdockq_8Acutoff': 'pDockQ'}
    if metric1 in metrics_renamed:
        plt.xlabel(metrics_renamed[metric1])
    else:
        plt.xlabel(metric1)

    if metric2 in metrics_renamed:
        plt.ylabel(metrics_renamed[metric2])
    else:
        plt.ylabel(metric2)

    # add correlation coefficient to plot
    plt.figtext(0.15, 0.85, '{} ρ = {}'.format(correlation_type.capitalize(), str(round(correlation, 2))))

    if output_filepath is not None:
        plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_capri_class_counts(results_df, dockq_col, output_dir, output_format='png'):
    # % of entries (or counts) in each CAPRI class (as defined by the thresholds in the DockQ github repo)
    if dockq_col == 'dockq':
        capri_col = 'CAPRI class'
    else:
        capri_col = '{} CAPRI class'.format(dockq_col.split('_')[:-1])
    df = copy.deepcopy(results_df).dropna(subset=dockq_col, axis=0)

    # create a col with CAPRI class based on DockQ thresholds
    df[capri_col] = df[dockq_col].apply(get_capri_class)

    palette=sns.color_palette('mako_r')
    plt.figure(figsize=(4, 4))
    sns.countplot(df, x=capri_col, order=['Incorrect', 'Acceptable', 'Medium', 'High'], palette=palette)
    plt.xlabel('CAPRI class')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, 'capri_countplot.{}'.format(output_format)), dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_comparison(df1, df2, column, metric, output_filepath, plot_type='violinplot', hue=None, color=None,
                    figsize=(4,4), fill=True):
    plt.figure(figsize=figsize)
    both_df = pd.concat([df1, df2], ignore_index=True).dropna(subset=metric, axis=0)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ', 'iptm': 'ipTM'}
    both_df.rename(mapper=metrics_renamed, axis=1, errors='ignore', inplace=True)
    if plot_type == 'boxplot':
        g = sns.boxplot(both_df, x=column, y=metrics_renamed[metric], hue=hue, color=color, fill=fill)
    elif plot_type == 'violinplot':
        g = sns.violinplot(both_df, x=column, y=metrics_renamed[metric], hue=hue, color=color, fill=fill,
                           width=0.6, cut=0)

    sns.despine()

    # add statistical annotations
    df1_data = df1[metric].dropna()
    df2_data = df2[metric].dropna()
    U, p = mannwhitneyu(df1_data,
                        df2_data,
                        alternative='two-sided')
    x1, x2 = 0, 1
    y = both_df[metrics_renamed[metric]].max() + 0.1
    h = 0.05
    plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c='black')
    if p < 0.001:
        sig_symbol = '***'
    elif p < 0.01:
        sig_symbol = '**'
    elif p < 0.05:
        sig_symbol = '*'
    else:
        sig_symbol = 'ns'
    print(f'p-value: {p}')
    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_comparison_virus_vs_bacteria(virus_results_df, bacteria_results_df, metric, output_filepath,
                                      plot_type='boxplot'):
    virus_results_df['Dataset'] = 'Virus' #'Virus-Mammalia Dimers'
    bacteria_results_df['Dataset'] = 'Bacteria' #'Bacteria-Mammalia Dimers'
    plot_comparison(df1=virus_results_df, df2=bacteria_results_df, column='Dataset',
                    metric=metric, output_filepath=output_filepath, plot_type=plot_type,
                    color=sns.color_palette('Paired')[0], figsize=(3.5, 3.5))


def plot_comparison_releasedate(results_df, metadata_df, metric, output_filepath, plot_type='boxplot'):
    # Analyze results by PDB release date (structures AF was trained on vs. structures released later)
    results_df['name'] = results_df['name'].str.replace('_fixedmodel', '')
    merged_df = results_df.merge(metadata_df, how='inner', left_on='name', right_on='assembly_id')
    merged_df['Release Date'] = merged_df['entry.rcsb_accession_info.initial_release_date'].apply(classify_release_date)
    before_df = merged_df[merged_df['Release Date'] == 'Before training set cutoff']
    after_df = merged_df[merged_df['Release Date'] == 'After training set cutoff']
    plot_comparison(df1=before_df, df2=after_df, column='Release Date',
                    metric=metric, output_filepath=output_filepath, plot_type=plot_type)


def plot_comparison_antibodies(results_df, metadata_df, metric, output_filepath, plot_type='boxplot'):
    results_df['name'] = results_df['name'].str.replace('_fixedmodel', '')
    merged_df = results_df.merge(metadata_df, how='inner', left_on='name', right_on='assembly_id')
    merged_df['Type of Entry'] = merged_df.apply(classify_antibody, axis=1)
    antibody_df = merged_df[merged_df['Type of Entry'] == 'Antibody']
    print(f'Entries with antibodies: {antibody_df.shape[0]}')
    notantibody_df = merged_df[merged_df['Type of Entry'] == 'Not antibody']
    print(f'Entries without antibodies: {notantibody_df.shape[0]}')
    plot_comparison(df1=antibody_df, df2=notantibody_df, column='Type of Entry',
                    metric=metric, output_filepath=output_filepath, plot_type=plot_type,
                    figsize=(3.5, 3.5), color=sns.color_palette('Paired')[0])


def plot_comparison_paired_test(df1, df2, column, metric, output_filepath, plot_type='violinplot',
                                figsize=(4,4)):
    plt.figure(figsize=figsize)
    df1 = df1.dropna(subset=metric, axis=0)
    df2 = df2.dropna(subset=metric, axis=0)
    pairs_in_both = list(set(df1['name'].to_list()).intersection(set(df2['name'].to_list())))
    df1_filtered = df1[df1['name'].isin(pairs_in_both)].sort_values(by=['name'])
    df2_filtered = df2[df2['name'].isin(pairs_in_both)].sort_values(by=['name'])
    both_df = pd.concat([df1, df2], ignore_index=True)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ', 'iptm': 'ipTM'}
    both_df.rename(mapper=metrics_renamed, axis=1, errors='ignore', inplace=True)

    if plot_type == 'boxplot':
        g = sns.boxplot(both_df, x=column, y=metrics_renamed[metric], hue=column)
    elif plot_type == 'violinplot':
        g = sns.violinplot(both_df, x=column, y=metrics_renamed[metric], hue=column, width=0.6, cut=0)

    sns.despine()

    # add statistical annotations
    df1_data = df1_filtered[metric]
    df2_data = df2_filtered[metric]
    diff = np.around(df1_data.values - df2_data.values, decimals=3)

    res = wilcoxon(diff)
    x1, x2 = 0, 1
    y = both_df[metrics_renamed[metric]].max() + 0.1
    h = 0.05
    plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c='black')
    if res.pvalue < 0.001:
        sig_symbol = '***'
    elif res.pvalue < 0.01:
        sig_symbol = '**'
    elif res.pvalue < 0.05:
        sig_symbol = '*'
    else:
        sig_symbol = 'ns'

    print(f'p-value: {res.pvalue}')

    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_comparison_input_type(pdbseqres_results_df, fullseq_results_df, metric, output_filepath, plot_type='boxplot'):
    pdbseqres_results_df['Dataset'] = 'PDB SEQRES'
    fullseq_results_df['Dataset'] = 'Full sequences'
    plot_comparison_paired_test(df1=pdbseqres_results_df, df2=fullseq_results_df,
                                column='Dataset', metric=metric, output_filepath=output_filepath,
                                plot_type=plot_type)


def plot_comparison_af_version(af2_results_df, af3_results_df, metric, output_filepath, plot_type='boxplot',
                               figsize=(3.5, 3.5)):
    af2_results_df['AlphaFold version'] = '2'
    af3_results_df['AlphaFold version'] = '3'
    plot_comparison_paired_test(df1=af2_results_df, df2=af3_results_df, column='AlphaFold version',
                                metric=metric, output_filepath=output_filepath, plot_type=plot_type,
                                figsize=figsize)
