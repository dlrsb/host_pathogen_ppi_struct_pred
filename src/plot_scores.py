import copy
import os
import datetime
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, mannwhitneyu, wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(font_scale=1.5)
# sns.set_style("white")
# SMALL_SIZE = 14
# MEDIUM_SIZE = 16
# BIGGER_SIZE = 18
# plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# plt.figure(figsize=(7.2, 7.2))
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


def plot_score_distribution(results_df, metric, output_dir, output_format='png', plot_type='hist_with_kde', bins=10):
    # plot types available: 'hist', 'kde', 'hist_with_kde
    data = results_df[metric].dropna()

    if plot_type == 'hist':
        sns.displot(x=data, bins=bins, height=4, aspect=1) # height 7.2
    elif plot_type == 'kde':
        sns.displot(x=data, kind='kde', height=4, aspect=1)
    elif plot_type == 'hist_with_kde':
        sns.displot(x=data, kde=True, bins=bins, height=4, aspect=1)

    # metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ'}
    #
    # if metric in metrics_renamed:
    #     plt.xlabel(metrics_renamed[metric])
    if metric == 'mmalign_tmscore':
        plt.xlabel('TM-score')
    elif metric == 'dockq':
        plt.xlabel('DockQ')
    plt.savefig(os.path.join(output_dir, '{}_{}.{}'.format(metric, plot_type, output_format)), dpi=300,
                bbox_inches='tight')
    plt.show()
    plt.close()


def plot_true_vs_predicted_scores(results_df, metric1, metric2, output_dir, output_format='png',
                                  correlation_type='pearson', figsize=(4, 4)):
    metrics_df = results_df[[metric1, metric2]].dropna(axis=0)
    metric1_data = metrics_df[metric2] # TODO: why do I change the order here? Fix this
    metric2_data = metrics_df[metric1]
    if correlation_type == 'pearson':
        correlation = pearsonr(metric1_data, metric2_data)[0]
    elif correlation_type == 'spearman':
        correlation = spearmanr(metric1_data, metric2_data)[0]

    # create scatterplot
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(data=results_df, x=metric1_data,
                    y=metric2_data, ax=ax, alpha=0.5)

    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ',
                       'jaccard_index': 'Jaccard Index',
                       'ialign8_itmscore': 'iAlign iTM-score (8 Å cutoff)'}
    if metric2 in metrics_renamed:
        print('here')
        plt.xlabel(metrics_renamed[metric2])
    else:
        plt.xlabel(metric2)

    if metric1 in metrics_renamed:
        plt.ylabel(metrics_renamed[metric1])
    else:
        plt.ylabel(metric1)

    # add correlation coefficient to plot

    plt.figtext(0.15, 0.85, '{} ρ = {}'.format(correlation_type.capitalize(), str(round(correlation, 2))))
    plt.savefig(os.path.join(output_dir, '{}_{}_scatterplot_{}.{}'.format(metric2, metric1,
                                                                          correlation_type, output_format)).replace("\\","/"),
                dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def corrfunc(x, y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = pearsonr(x, y) # TODO: see if I should change this to Spearman correlation
    ax = ax or plt.gca()
    ax.annotate(f'ρ = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)


def plot_scores_pairplot(results_df, metrics, output_dir, output_format='png', diag_kind='hist', corner=False):
    # diag_kind can be 'hist', 'kde' or 'auto
    # kind{‘scatter’, ‘kde’, ‘hist’, ‘reg’}
    # https://stackoverflow.com/questions/63416894/correlation-values-in-pairplot
    metrics_df = results_df[metrics].dropna(axis=0)
    # g = sns.pairplot(metrics_df, x_vars=x_metrics, y_vars=y_metrics, diag_kind=diag_kind, corner=corner)
    g = sns.pairplot(metrics_df, diag_kind=diag_kind, corner=corner)
    g.map_lower(corrfunc)
    # g.map_lower(sns.regplot)
    plt.savefig(os.path.join(output_dir, 'pairplot.{}'.format(output_format)))
    plt.show()
    plt.close()
    # TODO add stuff to change figure size/resolution, save as TIFF or SVG, etc.


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
    plt.savefig(os.path.join(output_dir, 'capri_countplot.{}'.format(output_format)), dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_comparison_virus_vs_bacteria(virus_results_df, bacteria_results_df, metric, output_dir, output_format='png',
                                      plot_type='boxplot'):
    virus_results_df['Dataset'] = 'Virus' #'Virus-Mammalia Dimers'
    bacteria_results_df['Dataset'] = 'Bacteria' #'Bacteria-Mammalia Dimers'
    both_df = pd.concat([virus_results_df, bacteria_results_df], ignore_index=True).dropna(subset=metric, axis=0)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ'}
    both_df.rename(mapper=metrics_renamed, axis=1, errors='ignore', inplace=True)

    if plot_type == 'boxplot':
        g = sns.boxplot(both_df, x='Dataset', y=metrics_renamed[metric])
    elif plot_type == 'violinplot':
        g = sns.violinplot(both_df, x='Dataset', y=metrics_renamed[metric], cut=0)

    # add statistical annotations
    virus_data = virus_results_df[metric].dropna()
    bacteria_data = bacteria_results_df[metric].dropna()
    U, p = mannwhitneyu(virus_data,
                        bacteria_data,
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
    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(os.path.join(output_dir, 'virus_vs_bacteria_{}_{}.{}'.format(metric, plot_type, output_format)))
    plt.show()
    plt.close()


def classify_release_date(release_date):
    date_str = release_date.split('T')[0]
    y, m, d = [int(x) for x in date_str.split('-')]
    date = datetime.date(y, m, d)
    cutoff_date = datetime.date(2021, 9, 30)
    if date > cutoff_date:
        return 'After training set cutoff'
    else:
        return 'Before training set cutoff'


def plot_comparison_releasedate(results_df, metadata_df, metric, output_dir, output_format='png', plot_type='boxplot'):
    # Analyze results by PDB release date (structures Alphafold-multimer v.2.3.1 was trained on vs.
    # structures released later)
    # TODO: confirm that merge is being done correctly
    # TODO: confirm that before and after groups are being separated correctly

    results_df['name'] = results_df['name'].str.replace('_fixedmodel', '')
    merged_df = results_df.merge(metadata_df, how='inner', left_on='name', right_on='assembly_id')
    merged_df['Release Date'] = merged_df['entry.rcsb_accession_info.initial_release_date'].apply(classify_release_date)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ'}
    merged_df.rename(metrics_renamed, axis=1, inplace=True)

    if plot_type == 'boxplot':
        g = sns.boxplot(merged_df, x='Release Date', y=metrics_renamed[metric])
    elif plot_type == 'violinplot':
        g = sns.violinplot(merged_df, x='Release Date', y=metrics_renamed[metric], cut=0)

    # add statistical annotations
    data_before = merged_df[merged_df['Release Date'] == 'Before training set cutoff'][metrics_renamed[metric]].dropna()
    data_after = merged_df[merged_df['Release Date'] == 'After training set cutoff'][metrics_renamed[metric]].dropna()

    U, p = mannwhitneyu(data_before,
                        data_after,
                        alternative='two-sided')
    x1, x2 = 0, 1
    y = merged_df[metrics_renamed[metric]].max() + 0.1
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
    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(
        os.path.join(output_dir, 'before_vs_after_training_cutoff_{}_{}.{}'.format(metric, plot_type, output_format)))
    plt.show()
    plt.close()


def plot_comparison_peptides_proteins(results_df, dataset_df, metric, output_dir, output_format='png',
                                      plot_type='boxplot'):
    # boxplots/violinplots side by side
    merged_df = results_df.merge(dataset_df, left_on='name', right_on='assembly_id')
    # Group by a categorical variable, referencing columns in a dataframe:
    # sns.boxplot(data=titanic, x="age", y="class")
    # sns.violinplot(data=df, x="age", y="class")
    # https://stackoverflow.com/questions/36578458/how-does-one-insert-statistical-annotations-stars-or-p-values
    # maybe concat both virus and bacteria datasets to analyze this


def classify_antibody(row):
    description_cols = ['polymer_entity_instances.1.polymer_entity.rcsb_polymer_entity.pdbx_description',
                        'polymer_entity_instances.2.polymer_entity.rcsb_polymer_entity.pdbx_description']
    antibody_nanobody_terms = ['nanobody', 'antibody', 'intrabody', 'vhh', 'nb8194', 'nb8193', 'mab 4e11', 'scfv',
                               'single chain variable fragment', 'single-chain fv', 'h11-h4', 'nm1230',
                               'nb_1A7', 'nb_1b11', 'single-domain antibody', 'single domain antibody',
                               'nm1226', 'nb112', 'nb1A2', 'nb2b4', 'n3113', 'nb-007', 'Ab08',
                               'nb70', 'vh ab6', 'immunoglobulin heavy chain variable region',
                               'nb113', 'ca1698', 'nb179', 'ca1697', 'variable fragment',
                               'cab-h7s', 'cab-f11n', 'nbfedf6', 'nbfedf7', 'variable region',
                               'ig gamma-3 chain c region', '5d', 'e3', '7f', 'nanoe6', 'nanof7', 'nanob12',
                               'jli-h11', 'jlk-g12', 'jli-g10', 'jle-e5', 'nb_msba 1', 'lam8']
    if any(term in row[description_cols[0]].lower() for term in antibody_nanobody_terms) or any(term in row[description_cols[1]].lower() for term in antibody_nanobody_terms):
        return 'Antibody'
    else:
        return 'Not antibody'


def plot_comparison_antibodies(results_df, metadata_df, metric, output_dir, output_format='png', plot_type='boxplot'):
    results_df['name'] = results_df['name'].str.replace('_fixedmodel', '')
    merged_df = results_df.merge(metadata_df, how='inner', left_on='name', right_on='assembly_id')
    merged_df['Type of Entry'] = merged_df.apply(classify_antibody, axis=1)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ'}
    merged_df.rename(metrics_renamed, axis=1, inplace=True)

    if plot_type == 'boxplot':
        g = sns.boxplot(merged_df, x='Type of Entry', y=metrics_renamed[metric])
    elif plot_type == 'violinplot':
        g = sns.violinplot(merged_df, x='Type of Entry', y=metrics_renamed[metric], cut=0)

    # add statistical annotations
    data_antibody = merged_df[merged_df['Type of Entry'] == 'Antibody'][metrics_renamed[metric]].dropna()
    data_notantibody = merged_df[merged_df['Type of Entry'] == 'Not antibody'][metrics_renamed[metric]].dropna()

    U, p = mannwhitneyu(data_antibody,
                        data_notantibody,
                        alternative='two-sided')
    x1, x2 = 0, 1
    y = merged_df[metrics_renamed[metric]].max() + 0.1
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
    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(
        os.path.join(output_dir, 'antibody_vs_nonantibody_{}_{}.{}'.format(metric, plot_type, output_format)))
    plt.show()
    plt.close()


def plot_comparison_antibodies_releasedate(results_df, metadata_df, metric, output_dir, output_format='png',
                                           plot_type='boxplot'):
    results_df['name'] = results_df['name'].str.replace('_fixedmodel', '')
    merged_df = results_df.merge(metadata_df, how='inner', left_on='name', right_on='assembly_id')
    merged_df['release_date'] = merged_df['entry.rcsb_accession_info.initial_release_date'].apply(classify_release_date)
    merged_df['type_of_entry'] = merged_df.apply(classify_antibody, axis=1)

    # Test homogeneity of variances

    # # Two-way ANOVA
    # model = ols('{} ~ C(release_date) + C(type_of_entry) + C(release_date):C(type_of_entry)'.format(metric),
    #             data=merged_df).fit()
    # # Test normality of residuals
    # shapiro_result = shapiro(model.resid)
    # print(shapiro_result)
    # jarque_bera_result = jarque_bera(model.resid)
    # print(jarque_bera_result)
    # # ANOVA
    # anova_results = sm.stats.anova_lm(model, typ=2)
    # #print(anova_results)
    # # Test normality of residuals
    # # Tukey's HSD

    cols_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ', 'release_date': 'Release Date',
                    'type_of_entry': 'Type of Entry'}
    merged_df.rename(cols_renamed, axis=1, inplace=True)
    fig, ax = plt.subplots(figsize=(11.7, 8.27))
    # Boxplot
    if plot_type == 'boxplot':
        g = sns.boxplot(merged_df, x='Type of Entry', y=cols_renamed[metric], hue='Release Date', ax=ax)
    elif plot_type == 'violinplot':
        g = sns.violinplot(merged_df, x='Type of Entry', y=cols_renamed[metric], hue='Release Date', cut=0, ax=ax)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.legend(loc='lower center')
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, 'antibody_and_releasedate_{}_{}.{}'.format(metric, plot_type, output_format)))
    plt.show()
    plt.close()
    # TODO: 2 factors


def plot_comparison_input_type(pdbseqres_results_df, fullseq_results_df, metric, output_dir, output_format='png', plot_type='boxplot'):
    pdbseqres_results_df['Dataset'] = 'PDB SEQRES'
    fullseq_results_df['Dataset'] = 'Full sequences'
    pdbseqres_results_df = pdbseqres_results_df.dropna(subset=metric, axis=0)
    fullseq_results_df = fullseq_results_df.dropna(subset=metric, axis=0)
    pairs_in_both = list(set(pdbseqres_results_df['name'].to_list()).intersection(set(fullseq_results_df['name'].to_list())))
    pdbseqres_results_df_filtered = pdbseqres_results_df[pdbseqres_results_df['name'].isin(pairs_in_both)].sort_values(by=['name'])
    fullseq_results_df_filtered = fullseq_results_df[fullseq_results_df['name'].isin(pairs_in_both)].sort_values(by=['name'])
    both_df = pd.concat([pdbseqres_results_df_filtered, fullseq_results_df_filtered], ignore_index=True)
    metrics_renamed = {'mmalign_tmscore': 'TM-score', 'dockq': 'DockQ'}
    both_df.rename(mapper=metrics_renamed, axis=1, errors='ignore', inplace=True)

    if plot_type == 'boxplot':
        g = sns.boxplot(both_df, x='Dataset', y=metrics_renamed[metric])
    elif plot_type == 'violinplot':
        g = sns.violinplot(both_df, x='Dataset', y=metrics_renamed[metric], cut=0)

    # add statistical annotations
    pdbseqres_data = pdbseqres_results_df_filtered[metric]
    fullseq_data = fullseq_results_df_filtered[metric]
    diff = np.around(pdbseqres_data.values - fullseq_data.values, decimals=3)

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
    plt.text((x1 + x2) * .5, y + h, sig_symbol, ha='center', va='bottom', color='black')
    plt.savefig(os.path.join(output_dir, 'pdbseqres_vs_fullseq_{}_{}.{}'.format(metric, plot_type, output_format)))
    plt.show()
    plt.close()


def plot_pathogens_high_iptm():
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=14)  # legend fontsize
    plt.rc('figure', titlesize=18)  # fontsize of the figure title
    other_pathogens_results_df = pd.read_csv(
        '../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    df_selected = other_pathogens_results_df[other_pathogens_results_df['iptm'] >= 0.5]

    def get_pathogen(x):
        if 'SARS-CoV1' in x['name']:
            pathogen = 'SARS-CoV-1'
        elif 'SARS-CoV2' in x['name']:
            pathogen = 'SARS-CoV-2'
        elif 'Mtb' in x['name']:
            pathogen = 'M. tuberculosis'
        else:
            pathogen = x['name'].split('-')[0]
        return pathogen

    df_selected['pathogen'] = df_selected.apply(get_pathogen, axis=1)
    plt.figure(figsize=(8, 8))
    sns.countplot(y=df_selected['pathogen'], orient='h', color='silver',
                  order = df_selected['pathogen'].value_counts().index)
    plt.ylabel('Pathogen')
    plt.xlabel('Number of protein pairs with > ipTM >= 0.5')
    plt.savefig('C:\\Users\\dlrba\\Desktop\\pathogens_high_iptm.png', bbox_inches='tight')


def plot_pdockq_success_rate():
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=14)  # legend fontsize
    plt.rc('figure', titlesize=18)  # fontsize of the figure title
    other_pathogens_results_df = pd.read_csv(
        '../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    print(other_pathogens_results_df.shape)
    print(other_pathogens_results_df[other_pathogens_results_df['iptm'] >= 0.5].shape)
    print(other_pathogens_results_df[other_pathogens_results_df['huintaf2_pdockq_8Acutoff'] >= 0.23].shape)
    plt.figure(figsize=(6, 6))
    plot_capri_class_counts(other_pathogens_results_df, dockq_col='huintaf2_pdockq_8Acutoff',
                            output_dir='C:\\Users\\dlrba\\Desktop')


def plot_num_virus_structures():
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=14)  # legend fontsize
    plt.rc('figure', titlesize=18)  # fontsize of the figure title
    other_pathogens_results_df = pd.read_csv(
        '../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    # TODO: include non-published datasets as well??
    #other_pathogens_results_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_all_onlybestmodel.csv')
    def get_pathogen(x):
        if 'SARS-CoV1' in x['name']:
            pathogen = 'SARS-CoV-1'
        elif 'SARS-CoV2' in x['name']:
            pathogen = 'SARS-CoV-2'
        elif 'Mtb' in x['name']:
            pathogen = 'M. tuberculosis'
        elif 'ZIKVfp' in x['name']:
            pathogen = 'ZIKV'
        else:
            pathogen = x['name'].split('-')[0]
        return pathogen

    other_pathogens_results_df['pathogen'] = other_pathogens_results_df.apply(get_pathogen, axis=1)
    df_viral = copy.deepcopy(other_pathogens_results_df[~other_pathogens_results_df['name'].str.contains('Chlamydia|Mtb', regex=True)])
    plt.figure(figsize=(5, 6))
    sns.countplot(y=df_viral['pathogen'], orient='h', color='silver',
                  order = df_viral['pathogen'].value_counts().index)
    plt.ylabel('Virus')
    plt.xlabel('Predicted structures')
    plt.savefig('C:\\Users\\dlrba\\Desktop\\virus_num_structures.png', dpi=300, bbox_inches='tight')


def plot_interface_comparison_scores(output_dir):
    df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/interface_analysis2/interface_comparison_results_all.csv')
    #print(df.shape[0])
    #print(df.columns)
    #print(df.isnull().sum(axis=0))
    # TODO: make this more general

    scores_pairs = [('jaccard_index', 'ialign_itmscore'),
                    ('jaccard_index', 'ialign_itmscore_rmsd'),
                    ('jaccard_index', 'ialign_isscore'),
                    ('jaccard_index', 'ialign_isscore_rmsd'),
                    ('jaccard_index', 'ialign8_itmscore'),
                    ('jaccard_index', 'ialign8_itmscore_rmsd'),
                    ('jaccard_index', 'ialign8_isscore'),
                    ('jaccard_index', 'ialign8_isscore_rmsd'),
                    ('jaccard_index', 'pcalign_pcscore'),
                    ('jaccard_index', 'foldseek_complexqtmscore'),
                    ('jaccard_index', 'foldseek_complexttmscore'),
                    #('jaccard_index', 'tmalign_TM-score_interface'),
                    #('jaccard_index', 'tmalign_RMSD_interface'),
                    ('ialign_itmscore', 'ialign_itmscore_rmsd'),
                    ('ialign_itmscore', 'pcalign_pcscore'),
                    ('ialign_itmscore', 'ialign_isscore'),
                    ('ialign_itmscore', 'foldseek_complexqtmscore'),
                    ('ialign_itmscore', 'foldseek_complexttmscore'),
                    ('ialign_isscore', 'ialign_isscore_rmsd'),
                    ('ialign_isscore', 'pcalign_pcscore'),
                    ('ialign_isscore', 'foldseek_complexqtmscore'),
                    ('ialign_isscore', 'foldseek_complexttmscore'),
                    ('ialign8_itmscore', 'ialign8_itmscore_rmsd'),
                    ('ialign8_itmscore', 'pcalign_pcscore'),
                    ('ialign8_itmscore', 'ialign8_isscore'),
                    ('ialign8_itmscore', 'foldseek_complexqtmscore'),
                    ('ialign8_itmscore', 'foldseek_complexttmscore'),
                    ('ialign8_isscore', 'ialign8_isscore_rmsd'),
                    ('ialign8_isscore', 'pcalign_pcscore'),
                    ('ialign8_isscore', 'foldseek_complexqtmscore'),
                    ('ialign8_isscore', 'foldseek_complexttmscore'),
                    ('pcalign_pcscore', 'foldseek_complexqtmscore'),
                    ('pcalign_pcscore', 'foldseek_complexttmscore')]
    for pair in scores_pairs:
        plot_true_vs_predicted_scores(df, pair[0], pair[1], correlation_type='spearman', figsize=(6, 6),
                                      output_dir=output_dir)


def plot_pathogens_max_jaccard():
    # TODO: make this more general
    plt.figure(figsize=(6, 6))
    df = pd.read_csv("C:\\Users\\dlrba\\Desktop\\interface_comparison_scores\\interface_comparison_results_all2.csv")
    cols = ['pathogen_interactor', 'bait', 'jaccard_index']
    df = df[cols]
    df_grouped = df.groupby(by=['pathogen_interactor', 'bait']).max()
    df_grouped.rename({'jaccard_index': 'Highest Jaccard Index'}, axis=1, inplace=True)
    sns.histplot(data=df_grouped, x='Highest Jaccard Index', binwidth=0.1, color='silver')
    plt.savefig('C:\\Users\\dlrba\\Desktop\\highest_jaccard_histplot.png')
    #df_grouped.to_csv('C:\\Users\\dlrba\\Desktop\\jaccard_index_grouped.csv')


def plot_num_pathogen_structures():
    plt.rc('font', size=14)  # controls default text sizes
    plt.rc('axes', titlesize=14)  # fontsize of the axes title
    plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels
    plt.rc('legend', fontsize=14)  # legend fontsize
    plt.rc('figure', titlesize=18)  # fontsize of the figure title
    other_pathogens_results_df = pd.read_csv(
        '../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_published_datasets_onlybestmodel.csv')
    # TODO: include non-published datasets as well??
    #other_pathogens_results_df = pd.read_csv('../../host_pathogen_ppi_struct_pred/results/krogan_lab_datasets/only_best_model/results_all_onlybestmodel.csv')
    def get_pathogen(x):
        if 'SARS-CoV1' in x['name']:
            pathogen = 'SARS-CoV'
        elif 'SARS-CoV2' in x['name']:
            pathogen = 'SARS-CoV-2'
        elif 'Mtb' in x['name']:
            pathogen = 'M. tuberculosis'
        elif 'ZIKVfp' in x['name']:
            pathogen = 'ZIKV'
        else:
            pathogen = x['name'].split('-')[0]
        return pathogen

    other_pathogens_results_df['pathogen'] = other_pathogens_results_df.apply(get_pathogen, axis=1)
    df_selected = other_pathogens_results_df[other_pathogens_results_df['iptm'] >= 0.5]
    all_val_counts = other_pathogens_results_df['pathogen'].value_counts()
    selected_val_counts = df_selected['pathogen'].value_counts()
    plt.figure(figsize=(5, 6))
    # sns.countplot(y=other_pathogens_results_df['pathogen'], orient='h', color='silver',
    #               order = other_pathogens_results_df['pathogen'].value_counts().index)
    sns.set_color_codes("pastel")
    sns.barplot(x=all_val_counts, y=all_val_counts.index,
                label="Total", color="b")

    sns.set_color_codes("muted")
    sns.barplot(x=selected_val_counts, y=selected_val_counts.index,
                label="iptm ≥ 0.5", color="b")

    plt.ylabel('Pathogen')
    plt.xlabel('Number of predicted structures')
    plt.legend(ncol=2, loc="lower right", frameon=True)
    plt.savefig('C:\\Users\\dbaptista\\Desktop\\pathogen_num_structures.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    plot_num_pathogen_structures()

