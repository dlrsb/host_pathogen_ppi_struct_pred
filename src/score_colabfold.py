import os
import glob
import json

from score_prediction import save_results
from scoring_metrics import pdockq_huintaf2_score


def get_scores_from_json(scores_filepath):
    scores_to_keep = ['ptm', 'iptm']
    with open(scores_filepath, 'r') as f:
        scores = json.load(f)
    selected_scores = {k: scores[k] for k in scores_to_keep}
    selected_scores['ranking_confidence'] = 0.8 * selected_scores['iptm'] + 0.2 * selected_scores['ptm']
    return selected_scores


def score_colabfold_structure(predictions_dir, interface_cutoff, output_filepath):
    # find best model:
    best_model_pth = glob.glob(os.path.join(predictions_dir, '*_unrelaxed_rank_001_*.pdb'))[0]
    print(best_model_pth)
    results = {'name': os.path.normpath(predictions_dir).split(os.sep)[-1],
               'ranked0_model_path': best_model_pth}
    best_model_scores_pth = glob.glob(os.path.join(predictions_dir, '*_scores_rank_001_*.json'))[0]
    af_results = get_scores_from_json(best_model_scores_pth)
    results.update(af_results)

    # calculate pDockQ as defined in the "Towards a structurally resolved human protein interaction network" paper,
    # with a user-defined cutoff value
    pdockq_huintaf2_results = pdockq_huintaf2_score(best_model_pth, cutoff=interface_cutoff)
    results.update(pdockq_huintaf2_results)

    # save results to file
    if output_filepath is None:
        output_filepath = os.path.join(predictions_dir, 'results.csv')
    save_results(results_dict=results, output_filepath=output_filepath)


def score_all(results_dir, output_filepath, interface_cutoff=8):
    prediction_dirs = [pred_dir for pred_dir in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, pred_dir))]
    for pred_dir in prediction_dirs:
        pred_dir_path = os.path.join(results_dir, pred_dir)
        try:
            score_colabfold_structure(pred_dir_path, interface_cutoff=interface_cutoff, output_filepath=output_filepath)
        except Exception as e:
            print('failed ' + pred_dir)
            print(e)


if __name__ == '__main__':
    score_all(results_dir='/cluster/project/beltrao/dbaptista/hpidb_af2_models/colabfold_predictions',
              output_filepath='/cluster/project/beltrao/dbaptista/hpidb_af2_models/hpidb_af2_results.csv',
              interface_cutoff=8)
