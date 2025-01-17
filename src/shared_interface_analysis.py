import copy
import os
import shutil
import ast
import subprocess
import itertools
import tempfile
import numpy as np
import scipy
import pandas as pd
from Bio.PDB import PDBParser
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scoring_metrics import rmsd_score_biopython_cealign
from rewrite_pdb import rewrite_pdb
from utils import load_structure, run_tmalign, create_temp_pdb, create_interface_temp_pdb, \
    select_residues_within_radius, run_pymol



def read_chain_len(fp_, chain_):
    # number of residues in a chain
    structure_ = PDBParser(QUIET=True).get_structure(fp_, fp_)
    return len([residue['CA'].coord for residue in structure_[0][chain_]]) #https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec203


def parse_resid(s):
    """ Return various list-of-residues representations as set-of-ints, e.g.:
    parse_resid('')
    parse_resid('1')
    parse_resid('1,2,3,4')
    parse_resid('{1,2,3,4}')
    parse_resid('[1,2,3,4]')
    """
    if isinstance(s, list) or isinstance(s, np.ndarray):
        return s
    if s == '' or s != s:
        return set()
    x = ast.literal_eval(s)
    if isinstance(x, int):
        return set([x])
    else:
        return set(map(int, x))


def jaccard_similarity(residues1, residues2):
    return len(set(residues1).intersection(set(residues2))) / len(set(residues1).union(set(residues2)))


def run_ialign(structure1_filepath, structure2_filepath, scoring_metric, normalization_method, distance_cutoff=4.5):
    # normalizing by average

    if scoring_metric == 'tm':
        score_name = 'TM-score'
    else:
        score_name = 'IS-score'
    score = np.nan
    rmsd = np.nan
    seqidentity = np.nan

    #output_path = os.path.join(workdir, 'outputs')
    #results_filepath = os.path.join(workdir, 'outputs', 'results.txt')
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_ialign.sh', '-s1',
             structure1_filepath, '-s2', structure2_filepath, '-n', normalization_method, '-m', scoring_metric,
             '--dc', str(distance_cutoff), '--minp', '5', '--mini', '3'], # previously '--minp', '10', '--mini', '8'
            check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        #print(lines)
        # f = open(results_filepath)
        # lines = [line.strip() for line in f.readlines()]
        # f.close()
        # shutil.rmtree(output_path)
        for line in lines:
            if score_name in line:
                score = line.split(' ')[2]
            elif 'RMSD' in line:
                split_line = line.split(' ')
                rmsd = split_line[3].replace(',', '')
                seqidentity = split_line[8]
    except subprocess.CalledProcessError as err:
        print(err)

    return score, rmsd, seqidentity


def run_intercomp(structure1_filepath, structure2_filepath, seed):
    intercomp_score = np.nan
    intercomp_seqscore = np.nan
    intercomp_pvalue = np.nan
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/tools/intercomp/InterfaceComparison', '-PDB', structure1_filepath,
             structure2_filepath, '-seed', seed], check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        for line in lines:
            fields = ['Best', 'sequence score', 'p-value']
            if all([x in line for x in fields]):
                split_line = line.split(' ')
                intercomp_score = split_line[1]
                intercomp_seqscore = split_line[4]
                intercomp_pvalue = split_line[12]
    except subprocess.CalledProcessError as err:
        print(err)

    return intercomp_score, intercomp_seqscore, intercomp_pvalue


def run_pcalign(structure1_filepath, structure2_filepath):
    pcalign_pcscore = np.nan

    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_pcalign.sh', '-s1',
             structure1_filepath, '-s2', structure2_filepath],
            check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        for line in lines:
            if 'PC-score' in line:
                pcalign_pcscore = line.split(' ')[2]
    except subprocess.CalledProcessError as err:
        print(err)

    return pcalign_pcscore


def run_foldseek_pairwise(structure1_filepath, structure2_filepath, structure1_interactor_chain, structure1_sharedprot_chain,
                          structure2_interactor_chain, structure2_sharedprot_chain, output_dir):

    output_path = os.path.join(output_dir, 'result')

    # initialize foldseek_results
    foldseek_results = {'foldseek_interactors_fident': np.nan,
                        'foldseek_interactors_evalue': np.nan,
                        'foldseek_interactors_bits': np.nan,
                        'foldseek_sharedprot_fident': np.nan,
                        'foldseek_sharedprot_evalue': np.nan,
                        'foldseek_sharedprot_bits': np.nan,
                        'foldseek_complexqtmscore': np.nan,
                        'foldseek_complexttmscore': np.nan}

    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_foldseek.sh', '-s1',
             structure1_filepath, '-s2', structure2_filepath, '--workdir', output_dir],
            check=True)
        df = pd.read_table(output_path, header=0)
        if df.shape[0] != 0:
            df['qchain'] = df['query'].str[-1]
            df['tchain'] = df['target'].str[-1]
            interactor_row = df[
                (df['qchain'] == structure1_interactor_chain) & (df['tchain'] == structure2_interactor_chain)]
            sharedprot_row = df[
                (df['qchain'] == structure1_sharedprot_chain) & (df['tchain'] == structure2_sharedprot_chain)]

            if (interactor_row.shape[0] != 0) and (sharedprot_row.shape[0] != 0):
                foldseek_results = {'foldseek_interactors_fident': interactor_row['fident'].item(),
                                    'foldseek_interactors_evalue': interactor_row['evalue'].item(),
                                    'foldseek_interactors_bits': interactor_row['bits'].item(),
                                    'foldseek_sharedprot_fident': sharedprot_row['fident'].item(),
                                    'foldseek_sharedprot_evalue': sharedprot_row['evalue'].item(),
                                    'foldseek_sharedprot_bits': sharedprot_row['bits'].item(),
                                    'foldseek_complexqtmscore': interactor_row['complexqtmscore'].item(),
                                    'foldseek_complexttmscore': interactor_row['complexttmscore'].item()}
                # complexqtmscore and complexttmscore are the same in both the interactor_row and sharedprot_row
            else:
                print(structure1_filepath)
                print(structure2_filepath)
                print('Chains matched by Foldseek are not the correct ones')
    except subprocess.CalledProcessError as err:
        print(err)

    return foldseek_results


def run_foldseek_all_vs_all(structures_dir, workdir):
    output_path = os.path.join(workdir, 'result')

    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_foldseek_all_vs_all.sh', '-d',
             structures_dir, '--workdir', workdir], check=True)
        foldseek_results = pd.read_table(output_path, header=0, sep='\t')
    except subprocess.CalledProcessError as err:
        print(err)
        foldseek_results = pd.DataFrame(columns=['query', 'target', 'fident', 'evalue', 'bits', 'complexqtmscore',
                                                 'complexttmscore', 'lddt', 'qtmscore', 'ttmscore', 'alntmscore',
                                                 'rmsd', 'prob'])

    return foldseek_results


def get_interactor(row, interactors_list):
    print(row)
    if row['interactor1'] in interactors_list:
        print(row['interactor1'])
        return row['interactor1']
    elif row['interactor2'] in interactors_list:
        print(row['interactor2'])
        return row['interactor2']


def run_ialign_list(pdb_list_filepath, scoring_metric, normalization_method, workdir_name, distance_cutoff=4.5):
    # normalizing by average
    workdir = '/cluster/project/beltrao/dbaptista/tools/ialign'

    if scoring_metric == 'tm':
        score_key = f'ialign{distance_cutoff}_itmscore'
    else:
        score_key = f'ialign{distance_cutoff}_isscore'

    results_dict = {'interaction_id1': [], 'interaction_id2': [], score_key: [], f'ialign{distance_cutoff}_rmsd': []}
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_ialign_list.sh', '-l',
             pdb_list_filepath, '-w', workdir_name, '-n', normalization_method, '-m', scoring_metric,
             '--dc', str(distance_cutoff), '--minp', '10', '--mini', '8'],
            check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        for i in range(len(lines)):
            if lines[i].startswith('>>>'):
                print(lines[i])
                line_split_i = lines[i].split(' ')
                results_dict['interaction_id1'].append(line_split_i[0][3:-2])
                results_dict['interaction_id2'].append(line_split_i[-1][0:-2])
                results_dict[score_key].append(lines[i + 5].split(' ')[2])
                split_line = lines[i + 8].split(' ')
                results_dict[f'ialign{distance_cutoff}_rmsd'].append(split_line[3].replace(',', ''))

    except subprocess.CalledProcessError as err:
        print(err)

    return results_dict


def run_ialign_pdb_vs_list(pdb_filepath, pdb_list_filepath, scoring_metric, normalization_method, workdir_name,
                           distance_cutoff=4.5):
    # normalizing by average

    if scoring_metric == 'tm':
        score_key = f'ialign{distance_cutoff}_itmscore'
    else:
        score_key = f'ialign{distance_cutoff}_isscore'

    results_dict = {'interaction_id1': [], 'interaction_id2': [], score_key: [], f'ialign{distance_cutoff}_rmsd': []}
    try:
        called_process = subprocess.run(
            ['/cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/scripts/run_ialign_pdb_vs_list.sh', '-p',
             pdb_filepath, '-l', pdb_list_filepath, '-w', workdir_name, '-n', normalization_method, '-m',
             scoring_metric, '--dc', str(distance_cutoff), '--minp', '10', '--mini', '8'],
            check=True, capture_output=True, text=True)
        lines = called_process.stdout.splitlines()
        for i in range(len(lines)):
            if lines[i].startswith('>>>'):
                print(lines[i])
                line_split_i = lines[i].split(' ')
                results_dict['interaction_id1'].append(line_split_i[0][3:-2])
                results_dict['interaction_id2'].append(line_split_i[-1][0:-2])
                results_dict[score_key].append(lines[i + 5].split(' ')[2])
                split_line = lines[i + 8].split(' ')
                results_dict[f'ialign{distance_cutoff}_rmsd'].append(split_line[3].replace(',', ''))

    except subprocess.CalledProcessError as err:
        print(err)

    return results_dict

class SharedInterfaceAnalysis():
    """
    def resid_updated_(r):
    res_A, res_B = calc_interface_residues(r.pdb, dist_threshold=5, plddt_threshold=0)
    if r.bait_chain == 'A':
        return res_A
    else:
        return res_B
    bait_.df_interactors['bait_ifresid'] = [ resid_updated_(r) for i, r in bait_.df_interactors.iterrows() ]
    """
    def __init__(self, df_interactors):
        self.uniprot_id_bait = df_interactors.head(1)['bait_id'].squeeze()
        print(self.uniprot_id_bait, 'bait uniprot id')
        self.df_interactors = df_interactors.copy()
        self.nresid_bait = read_chain_len(self.df_interactors.head(1).pdb.squeeze(), self.df_interactors.head(1).bait_chain.squeeze())
        print(self.nresid_bait, 'number of bait residues')
        self.df_interface_scores = None
        self.interactor_pairs = list(itertools.combinations(self.df_interactors['interactor_id'].to_list(), 2))
        self.df_matrix = None
        self.linkage_ = None

    def build_matrix(self):
        def apply_(r):
            return [resid in parse_resid(r['bait_ifresid']) for resid in range(1, self.nresid_bait + 1)] # this creates the binary vectors

        self.df_matrix = pd.DataFrame(self.df_interactors.apply(apply_, axis=1).to_list(), index=self.df_interactors['interactor_id'], columns=range(1, self.nresid_bait + 1))
        print(self.df_matrix.shape, 'shape matrix dimensions')
        print(sum(self.df_matrix), 'non-zero entries in residue matrix')

    def linkage(self, criterion='distance', t=.9):
        self.linkage_ = scipy.cluster.hierarchy.linkage(self.df_matrix, method='average', metric='jaccard')
        self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion=criterion, t=t)
        print(self.df_interactors['labels'].value_counts())

    def clustermap(self, fname=None):
        row_colors_ = [* self.df_interactors['labels'].map(lambda label: matplotlib.colormaps['tab10'].colors[label - 1]) ]
        # TODO: change colormap? in some cases, there may be more labels than the number of colors that 'tab10' provides
        if not (fname is None):
            plt.ioff()

        fontsize_pt = 14.0 # plt.rcParams['ytick.labelsize'] # medium = 10.0
        dpi = 72.27
        matrix_height_pt = fontsize_pt * self.df_matrix.shape[0]
        matrix_height_in = matrix_height_pt / dpi
        figure_height = matrix_height_in / (1 - 0.08)

        self.clustermap_ = sns.clustermap(self.df_matrix, row_linkage=self.linkage_, col_cluster=False, cbar_pos=None,
                                          figsize=(12, figure_height), row_colors=row_colors_, yticklabels=True) # added yticklabels
        # figsize was previously (12, 6)
        plt.title(self.uniprot_id_bait)
        if not (fname is None):
            #plt.savefig(fname, bbox_inches='tight', transparent=True)
            self.clustermap_.savefig(fname, bbox_inches='tight', transparent=True)
            plt.clf()
            plt.ion()

    def to_pymol(self, fname=None):
        """
            fname is:
                None => show pymol commands to execute locally
                file name => execute pymol on the cluster & store session as .pse
        """
        #df_interactors_top = self.df_interactors.sort_values('pdockq', ascending=False).groupby('labels').head(1) # Align to top interaction model
        df_interactors_top = self.df_interactors.sort_values(['labels', 'pdockq'], ascending=False)#.head(1) # Align to top interaction model
        cmd_ = []
        cmd_.append('delete all')
        cmd_.append('bg_color white')
        ref_ = df_interactors_top.head(1).squeeze() # Reference structure to align against
        for i, r in df_interactors_top.iterrows():
            if fname is None:
                fp_ = os.path.join('~/work-euler/', r.pdb.removeprefix('/cluster/work/beltrao/jjaenes/'))
            else:
                fp_ = r.pdb
                print (fp_)
            cmd_.append(f'load {fp_}, {r.interaction_id}')

            if r.interaction_id != ref_['interaction_id']:
                #cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')
                #chain_ = 'A' if r.interaction_id.startswith(r.bait_id) else 'B'
                #chain_ref_ = 'A' if ref_.interaction_id.startswith(ref_.bait_id) else 'B'

                cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')

            col_ = '0x' + matplotlib.colors.to_hex(matplotlib.colormaps['tab10'].colors[r.labels - 1])[1:]
            #interactor_chain = 'B' if r.bait_chain == 'A' else 'A'
            # if r.interaction_type == 'human-human':
            #     interactor_chain = 'B' if r.bait_chain == 'A' else 'A'
            # else:
            #     interactor_chain = 'C' if r.bait_chain == 'B' else 'B'
            interactor_chain = r.interactor_chain
            cmd_.append(f'color gray, {r.interaction_id} & chain {r.bait_chain}')
            cmd_.append(f'color {col_}, {r.interaction_id} & chain {interactor_chain}')

            #if r.interaction_id.startswith(ref_.bait_id):
            #    cmd_.append(f'color gray, {r.interaction_id} & chain A')
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain B')
            #else:
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain A')
            #    cmd_.append(f'color gray, {r.interaction_id} & chain B')

        for label, df_label in df_interactors_top.groupby('labels'):
            name = f'cluster{label}'
            print(name)
            members = ' '.join(df_label['interaction_id'])
            print(members)
            cmd_.append(f'group {name}, {members}')

        if fname is None:
            print('\n'.join(cmd_))
        else:
            cmd_.append(f'save {fname}')
            run_pymol(cmd_)

    def calc_jaccard_similarity(self, fname=None, only_compare_pathogen_against_human=False):
        jaccard_similarity = pd.Series(1 - scipy.spatial.distance.pdist(self.df_matrix.values, metric='jaccard'), name='jaccard_index',
                                          index=pd.MultiIndex.from_tuples(
                                              [(p1, p2) for i, p1 in enumerate(self.df_matrix.index.to_list()) for p2 in
                                               self.df_matrix.index[i + 1:]], names=['interactor1', 'interactor2']))
        jaccard_similarity_df = jaccard_similarity.reset_index()

        if only_compare_pathogen_against_human:
            pathogen_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'host-pathogen', 'interactor_id'].to_list()
            human_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'human-human', 'interactor_id'].to_list()
            pathogen_jaccard_df = copy.deepcopy(jaccard_similarity_df[jaccard_similarity_df['interactor1'].isin(pathogen_interactors) | jaccard_similarity_df['interactor2'].isin(pathogen_interactors)])
            pathogen_jaccard_df['pathogen_interactor'] = pathogen_jaccard_df.apply(get_interactor, interactors_list=pathogen_interactors, axis=1)
            pathogen_jaccard_df['human_interactor'] = pathogen_jaccard_df.apply(get_interactor, interactors_list=human_interactors, axis=1)
            pathogen_jaccard_df = pathogen_jaccard_df.drop(['interactor1', 'interactor2'], axis=1).dropna(axis=0, subset=['pathogen_interactor', 'human_interactor'])
            pathogen_jaccard_df = pathogen_jaccard_df[['pathogen_interactor', 'human_interactor', 'jaccard_index']]
            if fname is not None:
                pathogen_jaccard_df.to_csv(fname, index=False)
            return pathogen_jaccard_df
        else:
            if fname is not None:
                jaccard_similarity_df.to_csv(fname, index=False)
            return jaccard_similarity_df

    def calc_interface_comparison_scores(self):
        results_all = {'pathogen_interactor': [], 'human_interactor': [], 'ialign_itmscore': [], 'ialign_itmscore_rmsd': [],
                       'ialign_isscore': [], 'ialign_isscore_rmsd': [], 'ialign8_itmscore': [], 'ialign8_itmscore_rmsd': [],
                       'ialign8_isscore': [], 'ialign8_isscore_rmsd': [], 'pcalign_pcscore': []}
        pathogen_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'host-pathogen']
        human_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'human-human']
        for pathogen_interactor in pathogen_interactors['interactor_id'].to_list():
            struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'pdb'].item()
            for human_interactor in human_interactors['interactor_id'].to_list():
                struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == human_interactor, 'pdb'].item()

                results_all['pathogen_interactor'].append(pathogen_interactor)
                results_all['human_interactor'].append(human_interactor)

                print('calculating iAlign iTM-score (distance cutoff = 4.5)')
                itmscore, itmscore_rmsd, itmscore_seqidentity = run_ialign(struct1_filepath, struct2_filepath, scoring_metric='tm', normalization_method='ave')
                results_all['ialign_itmscore'].append(itmscore)
                results_all['ialign_itmscore_rmsd'].append(itmscore_rmsd)

                print('calculating iAlign IS-score (distance cutoff = 4.5)')
                isscore, isscore_rmsd, isscore_seqidentity = run_ialign(struct1_filepath, struct2_filepath, scoring_metric='is', normalization_method='ave')
                results_all['ialign_isscore'].append(isscore)
                results_all['ialign_isscore_rmsd'].append(isscore_rmsd)

                print('calculating iAlign iTM-score (distance cutoff = 8)')
                itmscore8, itmscore8_rmsd, itmscore8_seqidentity = run_ialign(struct1_filepath, struct2_filepath, scoring_metric='tm',
                                                                              normalization_method='ave',
                                                                              distance_cutoff=8)
                results_all['ialign8_itmscore'].append(itmscore8)
                results_all['ialign8_itmscore_rmsd'].append(itmscore8_rmsd)

                print('calculating iAlign IS-score (distance cutoff = 8)')
                isscore8, isscore8_rmsd, isscore8_seqidentity = run_ialign(struct1_filepath, struct2_filepath,
                                                                           scoring_metric='is',
                                                                           normalization_method='ave',
                                                                           distance_cutoff=8)
                results_all['ialign8_isscore'].append(isscore8)
                results_all['ialign8_isscore_rmsd'].append(isscore8_rmsd)

                # print('calculating intercomp score')
                # intercomp_score, intercomp_seqscore, intercomp_pvalue = run_intercomp(struct1_filepath, struct2_filepath, seed="12321")
                # results_all['intercomp_score'].append(intercomp_score)
                # results_all['intercomp_seqscore'].append(intercomp_seqscore)
                # results_all['intercomp_pvalue'].append(intercomp_pvalue)

                print('calculating pcalign pcscore')
                pcalign_pcscore = run_pcalign(struct1_filepath, struct2_filepath)
                results_all['pcalign_pcscore'].append(pcalign_pcscore)

        results_df = pd.DataFrame(results_all)
        if self.df_interface_scores is None:
            self.df_interface_scores = pd.DataFrame(results_all)
        else:
            scores_df = self.df_interface_scores.merge(results_df, how='inner', on=['pathogen_interactor', 'human_interactor'])
            self.df_interface_scores = scores_df

        return results_df

    def calc_tmscore(self, fname, workdir):
        results = {'pathogen_interactor': [], 'human_interactor': [], 'tmalign_tmscore': [], 'tmalign_rmsd': [],
                   'cealign_rmsd': []}
        pathogen_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'host-pathogen']
        human_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'human-human']
        #interactor_pairs = list(itertools.combinations(self.df_interactors['interactor_id'].to_list(), 2))
        with tempfile.TemporaryDirectory(dir=workdir) as tempdir:
            for pathogen_interactor in pathogen_interactors['interactor_id'].to_list():
                pathogen_interactor_info = pathogen_interactors[pathogen_interactors['interactor_id'] == pathogen_interactor]
                pathogen_interactor_pdb_path = create_temp_pdb(pathogen_interactor_info['pdb'].item(),
                                                               pathogen_interactor_info['interactor_chain'].item(), tempdir)

                for human_interactor in human_interactors['interactor_id'].to_list():
                    human_interactor_info = human_interactors[human_interactors['interactor_id'] == human_interactor]
                    human_interactor_pdb_path = create_temp_pdb(human_interactor_info['pdb'].item(),
                                                                human_interactor_info['interactor_chain'].item(), tempdir)

                    tmscore, tmalign_rmsd = run_tmalign(pathogen_interactor_pdb_path, human_interactor_pdb_path)

                    # CEalign
                    pathogen_interactor_structure = load_structure(pathogen_interactor_pdb_path)
                    human_interactor_structure = load_structure(human_interactor_pdb_path)
                    try:
                        cealign_rmsd = rmsd_score_biopython_cealign(human_interactor_structure, pathogen_interactor_structure)
                    except Exception as e:
                        print(e)
                        cealign_rmsd = np.nan

                    results['pathogen_interactor'].append(pathogen_interactor)
                    results['human_interactor'].append(human_interactor)
                    results['tmalign_tmscore'].append(tmscore)
                    results['tmalign_rmsd'].append(tmalign_rmsd)
                    results['cealign_rmsd'].append(cealign_rmsd)


            results_df = pd.DataFrame(data=results)
            results_df.to_csv(fname, index=False)

    def calc_tmscore_interface(self, radius, workdir):
        results_all = {'pathogen_interactor': [], 'human_interactor': [], 'tmalign_TM-score_interface': [],
                       'tmalign_RMSD_interface': []}
        pathogen_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'host-pathogen']
        human_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'human-human']

        with tempfile.TemporaryDirectory(dir=workdir) as tempdir:
            for pathogen_interactor in pathogen_interactors['interactor_id'].to_list():
                struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'pdb'].item()
                struct1_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'interactor_chain'].item()
                struct1_interface = ast.literal_eval(self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'interactor_ifresid'].item())

                residues_pathogen = select_residues_within_radius(pdb_filepath=struct1_filepath,
                                                         chain_id=struct1_chain_id,
                                                         residue_list=struct1_interface,
                                                         radius=radius)
                struct1_if_filepath = create_interface_temp_pdb(pdb_filepath=struct1_filepath,
                                                                chain_id=struct1_chain_id,
                                                                selected_residues=residues_pathogen,
                                                                output_dir=tempdir)
                for human_interactor in human_interactors['interactor_id'].to_list():
                    struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == human_interactor, 'pdb'].item()
                    struct2_chain_id = self.df_interactors.loc[
                        self.df_interactors['interactor_id'] == human_interactor, 'interactor_chain'].item()
                    struct2_interface = [int(x) for x in self.df_interactors.loc[
                        self.df_interactors['interactor_id'] == human_interactor, 'interactor_ifresid'].item().split(',')]

                    residues_human = select_residues_within_radius(pdb_filepath=struct2_filepath,
                                                             chain_id=struct2_chain_id,
                                                             residue_list=struct2_interface,
                                                             radius=radius)
                    struct2_if_filepath = create_interface_temp_pdb(pdb_filepath=struct2_filepath,
                                                                    chain_id=struct2_chain_id,
                                                                    selected_residues=residues_human,
                                                                    output_dir=tempdir)

                    tmscore, tmalign_rmsd = run_tmalign(struct1_if_filepath, struct2_if_filepath)
                    results_all['pathogen_interactor'].append(pathogen_interactor)
                    results_all['human_interactor'].append(human_interactor)
                    results_all['tmalign_TM-score_interface'].append(tmscore)
                    results_all['tmalign_RMSD_interface'].append(tmalign_rmsd)

        results_df = pd.DataFrame(results_all)
        if self.df_interface_scores is None:
            self.df_interface_scores = pd.DataFrame(results_all)
        else:
            scores_df = self.df_interface_scores.merge(results_df, how='inner', on=['pathogen_interactor', 'human_interactor'])
            self.df_interface_scores = scores_df

        return results_df

    def calc_foldseek_interface(self, workdir):
        results_all = {'pathogen_interactor': [], 'human_interactor': [],
                       'foldseek_interactors_fident': [], 'foldseek_interactors_evalue': [],
                       'foldseek_interactors_bits': [], 'foldseek_sharedprot_fident': [],
                       'foldseek_sharedprot_evalue': [], 'foldseek_sharedprot_bits': [],
                       'foldseek_complexqtmscore': [], 'foldseek_complexttmscore': []}
        pathogen_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'host-pathogen']
        human_interactors = self.df_interactors.loc[self.df_interactors['interaction_type'] == 'human-human']

        for pathogen_interactor in pathogen_interactors['interactor_id'].to_list():
            struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'pdb'].item()
            struct1_interactor_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'interactor_chain'].item()
            struct1_sharedprot_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == pathogen_interactor, 'bait_chain'].item()

            for human_interactor in human_interactors['interactor_id'].to_list():
                struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == human_interactor, 'pdb'].item()
                struct2_interactor_chain_id = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == human_interactor, 'interactor_chain'].item()
                struct2_sharedprot_chain_id = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == human_interactor, 'bait_chain'].item()

                foldseek_results = run_foldseek_pairwise(struct1_filepath, struct2_filepath,
                                                         structure1_interactor_chain=struct1_interactor_chain_id,
                                                         structure1_sharedprot_chain=struct1_sharedprot_chain_id,
                                                         structure2_interactor_chain=struct2_interactor_chain_id,
                                                         structure2_sharedprot_chain=struct2_sharedprot_chain_id,
                                                         output_dir=workdir)
                results_all['pathogen_interactor'].append(pathogen_interactor)
                results_all['human_interactor'].append(human_interactor)
                for key in foldseek_results:
                    results_all[key].append(foldseek_results[key])

        results_df = pd.DataFrame(results_all)
        if self.df_interface_scores is None:
            self.df_interface_scores = pd.DataFrame(results_all)
        else:
            scores_df = self.df_interface_scores.merge(results_df, how='inner', on=['pathogen_interactor', 'human_interactor'])
            self.df_interface_scores = scores_df

        return results_df

    def save_scores_to_file(self, fname):
        self.df_interface_scores.to_csv(fname, index=False)


class SharedInterfaceAnalysis2:
    """
    def resid_updated_(r):
    res_A, res_B = calc_interface_residues(r.pdb, dist_threshold=5, plddt_threshold=0)
    if r.bait_chain == 'A':
        return res_A
    else:
        return res_B
    bait_.df_interactors['bait_ifresid'] = [ resid_updated_(r) for i, r in bait_.df_interactors.iterrows() ]
    """
    def __init__(self, df_interactors):
        self.uniprot_id_bait = df_interactors.head(1)['bait_id'].squeeze()
        print(self.uniprot_id_bait, 'bait uniprot id')
        self.df_interactors = df_interactors.copy()
        self.nresid_bait = read_chain_len(self.df_interactors.head(1).pdb.squeeze(), self.df_interactors.head(1).bait_chain.squeeze())
        print(self.nresid_bait, 'number of bait residues')
        self.df_interface_scores = None
        self.interactor_pairs = list(itertools.combinations(self.df_interactors['interactor_id'].to_list(), 2))
        self.df_matrix = None
        self.linkage_ = None

    def build_matrix(self):
        def apply_(r):
            return [resid in parse_resid(r['bait_ifresid']) for resid in range(1, self.nresid_bait + 1)] # this creates the binary vectors

        self.df_matrix = pd.DataFrame(self.df_interactors.apply(apply_, axis=1).to_list(), index=self.df_interactors['interactor_id'], columns=range(1, self.nresid_bait + 1))
        print(self.df_matrix.shape, 'shape matrix dimensions')
        print(sum(self.df_matrix), 'non-zero entries in residue matrix')

    def linkage(self, criterion='distance', t=.9):
        self.linkage_ = scipy.cluster.hierarchy.linkage(self.df_matrix, method='average', metric='jaccard')
        self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion=criterion, t=t)
        print(self.df_interactors['labels'].value_counts())

    def clustermap(self, fname=None):
        row_colors_ = [* self.df_interactors['labels'].map(lambda label: matplotlib.colormaps['tab10'].colors[label - 1]) ]
        if not (fname is None):
            plt.ioff()

        fontsize_pt = 14.0 # plt.rcParams['ytick.labelsize'] # medium = 10.0
        dpi = 72.27
        matrix_height_pt = fontsize_pt * self.df_matrix.shape[0]
        matrix_height_in = matrix_height_pt / dpi
        figure_height = matrix_height_in / (1 - 0.08)

        self.clustermap_ = sns.clustermap(self.df_matrix, row_linkage=self.linkage_, col_cluster=False, cbar_pos=None,
                                          figsize=(12, figure_height), row_colors=row_colors_, yticklabels=True) # added yticklabels
        # figsize was previously (12, 6)
        plt.title(self.uniprot_id_bait)
        if not (fname is None):
            #plt.savefig(fname, bbox_inches='tight', transparent=True)
            self.clustermap_.savefig(fname, bbox_inches='tight', transparent=True)
            plt.clf()
            plt.ion()

    def to_pymol(self, fname=None):
        """
            fname is:
                None => show pymol commands to execute locally
                file name => execute pymol on the cluster & store session as .pse
        """
        #df_interactors_top = self.df_interactors.sort_values('pdockq', ascending=False).groupby('labels').head(1) # Align to top interaction model
        df_interactors_top = self.df_interactors.sort_values(['labels', 'pdockq'], ascending=False)#.head(1) # Align to top interaction model
        cmd_ = []
        cmd_.append('delete all')
        cmd_.append('bg_color white')
        ref_ = df_interactors_top.head(1).squeeze() # Reference structure to align against
        for i, r in df_interactors_top.iterrows():
            if fname is None:
                fp_ = os.path.join('~/work-euler/', r.pdb.removeprefix('/cluster/work/beltrao/jjaenes/'))
            else:
                fp_ = r.pdb
                print (fp_)
            cmd_.append(f'load {fp_}, {r.interaction_id}')

            if r.interaction_id != ref_['interaction_id']:
                #cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')
                #chain_ = 'A' if r.interaction_id.startswith(r.bait_id) else 'B'
                #chain_ref_ = 'A' if ref_.interaction_id.startswith(ref_.bait_id) else 'B'

                cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')

            col_ = '0x' + matplotlib.colors.to_hex(matplotlib.colormaps['tab10'].colors[r.labels - 1])[1:]
            #interactor_chain = 'B' if r.bait_chain == 'A' else 'A'
            # if r.interaction_type == 'human-human':
            #     interactor_chain = 'B' if r.bait_chain == 'A' else 'A'
            # else:
            #     interactor_chain = 'C' if r.bait_chain == 'B' else 'B'
            interactor_chain = r.interactor_chain
            cmd_.append(f'color gray, {r.interaction_id} & chain {r.bait_chain}')
            cmd_.append(f'color {col_}, {r.interaction_id} & chain {interactor_chain}')

            #if r.interaction_id.startswith(ref_.bait_id):
            #    cmd_.append(f'color gray, {r.interaction_id} & chain A')
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain B')
            #else:
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain A')
            #    cmd_.append(f'color gray, {r.interaction_id} & chain B')

        for label, df_label in df_interactors_top.groupby('labels'):
            name = f'cluster{label}'
            print(name)
            members = ' '.join(df_label['interaction_id'])
            print(members)
            cmd_.append(f'group {name}, {members}')

        if fname is None:
            print('\n'.join(cmd_))
        else:
            cmd_.append(f'save {fname}')
            run_pymol(cmd_)

    def calc_jaccard_similarity(self):
        jaccard_similarity = pd.Series(1 - scipy.spatial.distance.pdist(self.df_matrix.values, metric='jaccard'), name='jaccard_index',
                                          index=pd.MultiIndex.from_tuples(
                                              [(p1, p2) for i, p1 in enumerate(self.df_matrix.index.to_list()) for p2 in
                                               self.df_matrix.index[i + 1:]], names=['interactor1', 'interactor2']))
        jaccard_similarity_df = jaccard_similarity.reset_index()
        jaccard_similarity_df['bait'] = self.uniprot_id_bait

        self._add_to_scores_df(jaccard_similarity_df)

        return jaccard_similarity_df

    def calc_ialign_scores(self, scoring_metric, distance_cutoff):
        if scoring_metric == 'tm':
            scoring_metric_name = 'itmscore'
        elif scoring_metric == 'is':
            scoring_metric_name = 'isscore'

        results_all = {'interactor1': [], 'interactor2': [],
                       f'ialign{distance_cutoff}_{scoring_metric_name}': [],
                       f'ialign{distance_cutoff}_{scoring_metric_name}_rmsd': []}

        for interactor_pair in self.interactor_pairs:
            results_all['interactor1'].append(interactor_pair[0])
            results_all['interactor2'].append(interactor_pair[1])

            struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[0], 'pdb'].item()
            struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()

            print(f'calculating iAlign {scoring_metric_name} (distance cutoff = {distance_cutoff})')
            score, rmsd, seqidentity = run_ialign(struct1_filepath, struct2_filepath, scoring_metric=scoring_metric,
                                                  normalization_method='ave', distance_cutoff=distance_cutoff)
            results_all[f'ialign{distance_cutoff}_{scoring_metric_name}'].append(score)
            results_all[f'ialign{distance_cutoff}_{scoring_metric_name}_rmsd'].append(rmsd)

        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def calc_ialign_scores_all_vs_all(self, scoring_metric, distance_cutoff, ialign_dir):
        if scoring_metric == 'tm':
            scoring_metric_name = 'itmscore'
        elif scoring_metric == 'is':
            scoring_metric_name = 'isscore'

        # create pdb_list file
        pdb_filepaths = self.df_interactors['pdb'].to_list()
        interaction_ids = self.df_interactors['interaction_id'].to_list()
        new_pdb_filepaths = []
        for filepath, interaction_id in zip(pdb_filepaths, interaction_ids):
            new_filepath = os.path.join(ialign_dir, f'{interaction_id}.pdb')
            shutil.copyfile(filepath, new_filepath)
            fixed_filepath = os.path.join(ialign_dir, f'{interaction_id}_fixed.pdb')
            rewrite_pdb(new_filepath, fixed_filepath)
            new_pdb_filepaths.append(fixed_filepath)

        pdb_files_str = '\n'.join(new_pdb_filepaths)
        pdb_files_list_path = os.path.join(ialign_dir, f'{self.uniprot_id_bait}_pdb_files.lst')
        with open(pdb_files_list_path, 'w') as f:
            f.write(pdb_files_str)

        ialign_results = run_ialign_list(pdb_list_filepath=pdb_files_list_path, scoring_metric=scoring_metric,
                                         normalization_method='ave', distance_cutoff=distance_cutoff,
                                         workdir_name=self.uniprot_id_bait)

        def get_interactor_id(x):
            x_split = x.split('_')
            for id in x_split:
                if id != self.uniprot_id_bait:
                    return id

        results_df = pd.DataFrame(ialign_results)
        results_df['interactor1'] = results_df['interaction_id1'].apply(get_interactor_id)
        results_df['interactor2'] = results_df['interaction_id2'].apply(get_interactor_id)
        results_df.drop(['interaction_id1', 'interaction_id2'], axis=1, inplace=True)
        self._add_to_scores_df(results_df)

        for filepath in new_pdb_filepaths:
            if os.path.exists(filepath):
                os.remove(filepath)

        return results_df

    def calc_intercomp_score(self):
        results_all = {'interactor1': [], 'interactor2': [], 'intercomp_score': [],
                       'intercomp_seqscore': [], 'intercomp_pvalue': []}

        for interactor_pair in self.interactor_pairs:
            results_all['interactor1'].append(interactor_pair[0])
            results_all['interactor2'].append(interactor_pair[1])

            struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[0], 'pdb'].item()
            struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()

            intercomp_score, intercomp_seqscore, intercomp_pvalue = run_intercomp(struct1_filepath, struct2_filepath, seed="12321")
            results_all['intercomp_score'].append(intercomp_score)
            results_all['intercomp_seqscore'].append(intercomp_seqscore)
            results_all['intercomp_pvalue'].append(intercomp_pvalue)

        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def calc_pcalign_score(self):
        results_all = {'interactor1': [], 'interactor2': [], 'pcalign_pcscore': []}

        for interactor_pair in self.interactor_pairs:
            results_all['interactor1'].append(interactor_pair[0])
            results_all['interactor2'].append(interactor_pair[1])

            struct1_filepath = self.df_interactors.loc[
                self.df_interactors['interactor_id'] == interactor_pair[0], 'pdb'].item()
            struct2_filepath = self.df_interactors.loc[
                self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()

            pcalign_pcscore = run_pcalign(struct1_filepath, struct2_filepath)
            results_all['pcalign_pcscore'].append(pcalign_pcscore)

        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def calc_tmscore(self, workdir, interface_only=True, radius=5):
        if interface_only:
            results_all = {'interactor1': [], 'interactor2': [], 'tmalign_tmscore_interface': [],
                           'tmalign_rmsd_interface': []}
        else:
            results_all = {'interactor1': [], 'interactor2': [], 'tmalign_tmscore': [], 'tmalign_rmsd': []}
        with tempfile.TemporaryDirectory(dir=workdir) as tempdir:
            for interactor_pair in self.interactor_pairs:
                struct1_filepath = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == interactor_pair[0], 'pdb'].item()
                struct1_interactor_chain = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == interactor_pair[0], 'interactor_chain'].item()
                struct2_filepath = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()
                struct2_interactor_chain = self.df_interactors.loc[
                    self.df_interactors['interactor_id'] == interactor_pair[1], 'interactor_chain'].item()

                if interface_only:
                    struct1_interface = list(parse_resid(self.df_interactors.loc[
                                                             self.df_interactors['interactor_id'] == interactor_pair[
                                                                 0], 'interactor_ifresid'].item()))
                    interactor1_selected_residues = select_residues_within_radius(pdb_filepath=struct1_filepath,
                                                                                  chain_id=struct1_interactor_chain,
                                                                                  residue_list=struct1_interface,
                                                                                  radius=radius)
                    struct1_tmalign_input_file = create_interface_temp_pdb(pdb_filepath=struct1_filepath,
                                                                    chain_id=struct1_interactor_chain,
                                                                    selected_residues=interactor1_selected_residues,
                                                                    output_dir=tempdir)
                    struct2_filepath = self.df_interactors.loc[
                        self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()
                    struct2_interactor_chain = self.df_interactors.loc[
                        self.df_interactors['interactor_id'] == interactor_pair[1], 'interactor_chain'].item()
                    struct2_interface = list(parse_resid(self.df_interactors.loc[
                                                             self.df_interactors['interactor_id'] == interactor_pair[
                                                                 1], 'interactor_ifresid'].item()))
                    interactor2_selected_residues = select_residues_within_radius(pdb_filepath=struct2_filepath,
                                                                                  chain_id=struct2_interactor_chain,
                                                                                  residue_list=struct2_interface,
                                                                                  radius=radius)
                    struct2_tmalign_input_file = create_interface_temp_pdb(pdb_filepath=struct2_filepath,
                                                                    chain_id=struct2_interactor_chain,
                                                                    selected_residues=interactor2_selected_residues,
                                                                    output_dir=tempdir)
                    suffix = '_interface'
                else:
                    struct1_tmalign_input_file = create_temp_pdb(struct1_filepath, struct1_interactor_chain, tempdir)
                    struct2_tmalign_input_file = create_temp_pdb(struct2_filepath, struct2_interactor_chain, tempdir)
                    suffix = ''

                tmscore, tmalign_rmsd = run_tmalign(struct1_tmalign_input_file, struct2_tmalign_input_file)
                # TODO: change the scores i get from tmalign. The way it is now I'm using TM2, which considers struct2
                #  as the reference structure

                results_all['interactor1'].append(interactor_pair[0])
                results_all['interactor2'].append(interactor_pair[1])
                results_all[f'tmalign_tmscore{suffix}'].append(tmscore)
                results_all[f'tmalign_rmsd{suffix}'].append(tmalign_rmsd)


        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def calc_foldseek_interface(self, workdir):
        results_all = {'interactor1': [], 'interactor2': [],
                       'foldseek_interactors_fident': [], 'foldseek_interactors_evalue': [],
                       'foldseek_interactors_bits': [], 'foldseek_sharedprot_fident': [],
                       'foldseek_sharedprot_evalue': [], 'foldseek_sharedprot_bits': [],
                       'foldseek_complexqtmscore': [], 'foldseek_complexttmscore': []}

        for interactor_pair in self.interactor_pairs:
            results_all['interactor1'].append(interactor_pair[0])
            results_all['interactor2'].append(interactor_pair[1])

            struct1_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[0], 'pdb'].item()
            struct1_interactor_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[0], 'interactor_chain'].item()
            struct1_sharedprot_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[0], 'bait_chain'].item()
            struct2_filepath = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[1], 'pdb'].item()
            struct2_interactor_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[1], 'interactor_chain'].item()
            struct2_sharedprot_chain_id = self.df_interactors.loc[self.df_interactors['interactor_id'] == interactor_pair[1], 'bait_chain'].item()

            foldseek_results = run_foldseek_pairwise(struct1_filepath, struct2_filepath,
                                                     structure1_interactor_chain=struct1_interactor_chain_id,
                                                     structure1_sharedprot_chain=struct1_sharedprot_chain_id,
                                                     structure2_interactor_chain=struct2_interactor_chain_id,
                                                     structure2_sharedprot_chain=struct2_sharedprot_chain_id,
                                                     output_dir=workdir)

            for key in foldseek_results:
                results_all[key].append(foldseek_results[key])

        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def calc_foldseek_interface_all_vs_all(self, workdir):
        results_all = {'interactor1': [], 'interactor2': [],
                       'foldseek_interactors_fident': [], 'foldseek_interactors_evalue': [],
                       'foldseek_interactors_bits': [], 'foldseek_bait_fident': [],
                       'foldseek_bait_evalue': [], 'foldseek_bait_bits': [],
                       'foldseek_interactors_complexqtmscore': [], 'foldseek_interactors_complexttmscore': [],
                       'foldseek_bait_complexqtmscore': [], 'foldseek_bait_complexttmscore': []}

        pdb_filepaths = self.df_interactors['pdb'].to_list()
        interaction_ids = self.df_interactors['interaction_id'].to_list()
        pdb_dir = os.path.join(workdir, 'pdbfiles')
        if os.path.exists(pdb_dir):
            shutil.rmtree(pdb_dir)
        os.mkdir(pdb_dir)
        for filepath, interaction_id in zip(pdb_filepaths, interaction_ids):
            shutil.copyfile(filepath, os.path.join(pdb_dir, f'{interaction_id}.pdb'))

        foldseek_df = run_foldseek_all_vs_all(structures_dir=os.path.join(workdir, 'pdbfiles'), workdir=workdir)

        def get_info(row):
            return row['interaction_id'].item(), row['bait_chain'].item(), row['interactor_chain'].item()

        for pair in self.interactor_pairs:
            interaction1_id, struct1_bait_chain_id, struct1_interactor_chain_id = get_info(self.df_interactors[self.df_interactors['interactor_id'] == pair[0]])
            interaction2_id, struct2_bait_chain_id, struct2_interactor_chain_id = get_info(self.df_interactors[self.df_interactors['interactor_id'] == pair[1]])
            struct1_interactor_foldseek_id = f'{interaction1_id}.pdb_{struct1_interactor_chain_id}'
            struct2_interactor_foldseek_id = f'{interaction2_id}.pdb_{struct2_interactor_chain_id}'
            struct1_bait_foldseek_id = f'{interaction1_id}.pdb_{struct1_bait_chain_id}'
            struct2_bait_foldseek_id = f'{interaction2_id}.pdb_{struct2_bait_chain_id}'
            #res = foldseek_df[((foldseek_df['query'] == struct1_interactor_foldseek_id) & (foldseek_df['target'] == struct2_interactor_foldseek_id)) | ((foldseek_df['target'] == struct1_interactor_foldseek_id) & (foldseek_df['query'] == struct2_interactor_foldseek_id))]
            interactors_result = foldseek_df[(foldseek_df['query'] == struct1_interactor_foldseek_id) & (foldseek_df['target'] == struct2_interactor_foldseek_id)]
            bait_result = foldseek_df[(foldseek_df['query'] == struct1_bait_foldseek_id) & (foldseek_df['target'] == struct2_bait_foldseek_id)]

            foldseek_cols_to_keep = ['fident', 'evalue', 'bits', 'complexqtmscore', 'complexttmscore']
            results_all['interactor1'].append(pair[0])
            results_all['interactor2'].append(pair[1])

            if interactors_result.shape[0] == 0:
                interactors_result = pd.DataFrame(columns=foldseek_df.columns)
                interactors_result.loc[0, :] = np.nan
                interactors_result['query'] = struct1_interactor_foldseek_id
                interactors_result['target'] = struct2_interactor_foldseek_id

            if bait_result.shape[0] == 0:
                bait_result = pd.DataFrame(columns=foldseek_df.columns)
                bait_result.loc[0, :] = np.nan
                bait_result['query'] = struct1_bait_foldseek_id
                bait_result['target'] = struct2_bait_foldseek_id

            for col in foldseek_cols_to_keep:
                results_all[f'foldseek_interactors_{col}'].append(interactors_result[col].item())
                results_all[f'foldseek_bait_{col}'].append(bait_result[col].item())

        print(results_all)

        results_df = pd.DataFrame(results_all)
        self._add_to_scores_df(results_df)

        return results_df

    def _add_to_scores_df(self, results):
        results['bait'] = self.uniprot_id_bait
        if self.df_interface_scores is None:
            self.df_interface_scores = results
        else:
            scores_df = self.df_interface_scores.merge(results, how='left', on=['bait', 'interactor1', 'interactor2'])
            self.df_interface_scores = scores_df

    def save_scores_to_file(self, filepath):
        self.df_interface_scores.to_csv(filepath, index=False)