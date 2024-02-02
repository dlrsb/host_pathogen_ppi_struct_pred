import os
import glob
import json
import argparse
import pandas as pd
from rename_protein_chains import rename_protein_chains2
from msa_seqs import get_num_msa_hits
from select_subsequence_pdb import select_residues_pdb
from scoring_metrics import *


def save_results(results_dict, output_filepath):
    results_df = pd.DataFrame(data=results_dict, index=[0])
    # results_df.to_csv(output_filepath, index=False, mode='a', header=not os.path.exists(output_filepath)) # This way of saving seems to concatenate without taking column names into account. For example, when DockQ is not well calculated, the remaining columns are shifted
    if os.path.exists(output_filepath):
        saved_df = pd.read_csv(output_filepath)
        concat_df = pd.concat([saved_df, results_df], axis=0, ignore_index=True, sort=False)
        concat_df.to_csv(output_filepath, index=False)
    else:
        results_df.to_csv(output_filepath, index=False)


def add_prefix(results_dict, prefix):
    return {prefix + k: v for k, v in results_dict.items()}


def score_structure(predictions_dir, native_structure_filepath=None, input_fasta_filepath=None, rename_chains=False,
                    only_best_model=False, select_residues_before_dockq=False, interface_cutoff=10,
                    output_filepath=None):
    results = {'name': os.path.normpath(predictions_dir).split(os.sep)[-1]}

    # get all models in predictions_dir or only the best model
    models_to_score = []

    with open(os.path.join(predictions_dir, 'ranking_debug.json'), 'r') as f:
        ranking_debug = json.load(f)
        results['ranked0_model'] = ranking_debug['order'][0]

    # if os.path.exists(os.path.join(predictions_dir, 'relaxed_{}.pdb'.format(ranking_debug['order'][0]))):
    #     pdbs_to_use = 'relaxed'
    # else:
    #     pdbs_to_use = 'unrelaxed'
    pdbs_to_use = 'unrelaxed'

    if only_best_model:
        models_to_score.append('{}_{}.pdb'.format(pdbs_to_use, ranking_debug['order'][0]))
    else:
        pdb_paths = glob.glob(os.path.join(predictions_dir, pdbs_to_use + '_model_*[1-5]_multimer_v3_pred_*[0-5].pdb'))

        models_to_score.extend([os.path.split(pth)[-1] for pth in pdb_paths])

    for model in models_to_score:
        if not only_best_model:
            split_name = model.split('_')
            prefix = '{}_'.format('_'.join(split_name[1:3] + [split_name[5]] + [split_name[6].split('.')[0]]))
        else:
            prefix = ''

        # rename_chains
        if rename_chains:
            rename_protein_chains2(os.path.join(predictions_dir, model), input_fasta_filepath)
            model_filepath = os.path.join(predictions_dir, '{}_chains_renamed.pdb'.format(model.split('.')[0]))
        else:
            model_filepath = os.path.join(predictions_dir, model)

        # calculate scores for each structure (if native_structure_filepath exists, calculate TM-score and DockQ,
        # otherwise just calculate pDockQ)

        # get scores (ptm & iptm) provided by Alphafold-multimer
        af_scores_file = 'result_{}.pkl'.format('_'.join(model.split('.')[0].split('_')[1:]))
        if not os.path.exists(os.path.join(predictions_dir, af_scores_file)):
            af_scores_file = 'result_{}_reduced.pkl'.format('_'.join(model.split('.')[0].split('_')[1:]))
        af_results = get_alphafold_results_scores(os.path.join(predictions_dir, af_scores_file))
        results.update(add_prefix(af_results, prefix))

        # calculate pDockQ as defined in the Folddock paper, with the original 8A cutoff value:
        pdockq_folddock_results_8A = pdockq_folddock_score(model_filepath, cutoff=8)
        results.update(add_prefix(pdockq_folddock_results_8A, prefix))
        # # calculate pDockQ as defined in the Folddock paper, with the user-defined cutoff value:
        # pdockq_folddock_results = pdockq_folddock_score(model_filepath, cutoff=interface_cutoff)
        # results.update(add_prefix(pdockq_folddock_results, prefix))

        # calculate pDockQ as defined in the "Towards a structurally resolved human protein interaction network" paper,
        # with the original 10 A cutoff value:
        pdockq_huintaf2_results_10A = pdockq_huintaf2_score(model_filepath, cutoff=10)
        results.update(add_prefix(pdockq_huintaf2_results_10A, prefix))

        # calculate pDockQ as defined in the "Towards a structurally resolved human protein interaction network" paper,
        # with a user-defined cutoff value
        pdockq_huintaf2_results = pdockq_huintaf2_score(model_filepath, cutoff=interface_cutoff)
        results.update(add_prefix(pdockq_huintaf2_results, prefix))

        # calculate pDockQ2 as defined in "Evaluation of AlphaFold-Multimer prediction on multi-chain protein complexes"
        pdockq2_results_8 = pdockq2_score(predicted_structure_filepath=model_filepath,
                                          results_pickle_filepath=os.path.join(predictions_dir, af_scores_file),
                                          cutoff=8)
        results.update(add_prefix(pdockq2_results_8, prefix))

        # # calculate pDockQ2 using a user-defined cutoff value
        # pdockq2_results = pdockq2_score(predicted_structure_filepath=model_filepath,
        #                                 results_pickle_filepath=os.path.join(predictions_dir, af_scores_file),
        #                                 cutoff=interface_cutoff)
        # results.update(add_prefix(pdockq2_results, prefix))

        if native_structure_filepath is not None:

            tmscore_results = tm_score(native_structure_filepath=native_structure_filepath,
                                       predicted_structure_filepath=model_filepath)
            results.update(add_prefix(tmscore_results, prefix))
            if select_residues_before_dockq:
                select_residues_pdb(model_pdb_filepath=model_filepath,
                                    model_input_fasta=input_fasta_filepath,
                                    native_pdb_filepath=native_structure_filepath)
                new_model_filepath = '{}_selected.pdb'.format(os.path.splitext(model_filepath)[0])
                dockq_results = dockq_score(native_structure_filepath=native_structure_filepath,
                                            predicted_structure_filepath=new_model_filepath)
            else:
                dockq_results = dockq_score(native_structure_filepath=native_structure_filepath,
                                            predicted_structure_filepath=model_filepath)
            results.update(add_prefix(dockq_results, prefix))

        # calculate interface size for each structure
        # interface_results = get_interface_size(model_filepath, contacts_threshold=interface_cutoff)
        # results.update(add_prefix(interface_results, prefix))

    # get # of seqs in MSAs for each structure
    num_msas = get_num_msa_hits(predictions_dir)
    results.update(num_msas)

    # save results to file
    if output_filepath is None:
        output_filepath = os.path.join(predictions_dir, 'results.csv')
    save_results(results_dict=results, output_filepath=output_filepath)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Calculate scores, interface statistics and number of seqs in MSAs for AlphaFold2 predicted structures')
    arg_parser.add_argument('predictions_dir', help='Folder containing Alphafold-multimer outputs')
    arg_parser.add_argument(
        '-n',
        dest='native_structure_filepath',
        default=None,
        help='Path to the PDB file containing the native (experimentally-determined) structure',
    )
    arg_parser.add_argument(
        '-f',
        dest='input_fasta_filepath',
        default=None,
        help='Path to the input FASTA file (used when --rename-chains is true)',
    )
    arg_parser.add_argument(
        '-o',
        dest='output_filepath',
        default=None,
        help='Path to output file',
    )
    arg_parser.add_argument(
        '--interface-cutoff',
        dest='interface_cutoff',
        default=10,
        type=int,
        help='Distance cutoff value to define the interface',
    )
    arg_parser.add_argument(
        '--rename-chains',
        dest='rename_chains',
        action='store_true',
        help='Rename the protein chains in the predicted structure',
    )
    arg_parser.add_argument(
        '--only-best-model',
        dest='only_best_model',
        action='store_true',
        help='Only calculate results for the highest ranked model',
    )
    arg_parser.add_argument(
        '--select-residues-before-dockq',
        dest='select_residues_before_dockq',
        action='store_true',
        help='Only calculate DockQ based on the residues that also exist in the native file',
    )
    args = arg_parser.parse_args()
    score_structure(**vars(args))
