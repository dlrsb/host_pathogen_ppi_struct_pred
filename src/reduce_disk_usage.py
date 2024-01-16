import os
import glob
import shutil
import argparse
import pickle

from msa_seqs import get_num_msa_hits


def reduce_results_files(predictions_dir):
    results_pkl_paths = glob.glob(os.path.join(predictions_dir, 'result_model_*[1-5]_multimer_v3_pred_*[0-5].pkl'))
    scores_to_keep = ['num_recycles', 'predicted_aligned_error', 'plddt', 'ptm', 'iptm', 'ranking_confidence']
    # scores_to_keep = ['num_recycles', 'ptm', 'iptm', 'ranking_confidence']

    for pth in results_pkl_paths:
        with open(pth, 'rb') as f:
            af_results = pickle.load(f)
        selected_scores = {k: v for k, v in af_results.items() if k in scores_to_keep}
        new_pth = '{}_reduced{}'.format(os.path.splitext(pth)[0], os.path.splitext(pth)[1])
        with open(new_pth, 'wb') as f:
            pickle.dump(selected_scores, f)
        os.remove(pth)


def summarize_msa_hits(predictions_dir):
    msa_results = get_num_msa_hits(predictions_dir)
    with open(os.path.join(predictions_dir, 'msas', 'msa_hits.pkl'), 'wb') as f:
        pickle.dump(msa_results, f)

    for x in os.listdir(os.path.join(predictions_dir, 'msas')):
        pth = os.path.join(predictions_dir, 'msas', x)
        if os.path.isdir(pth):
            shutil.rmtree(pth)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Clean up AF2 results folder (only keeps the results.pkl file for the highest ranked model')
    arg_parser.add_argument('predictions_dir', help='Folder containing Alphafold-multimer outputs')
    arg_parser.add_argument(
        '--reduce-results',
        dest='reduce_results',
        action='store_true',
        help='Only keep selected scores in results pkl files',
    )
    arg_parser.add_argument(
        '--remove-msa-hits',
        dest='remove_msa_hits',
        action='store_true',
        help='Remove MSA hits',
    )
    args = arg_parser.parse_args()
    if args.reduce_results:
        reduce_results_files(args.predictions_dir)

    if args.remove_msa_hits:
        summarize_msa_hits(args.predictions_dir)


