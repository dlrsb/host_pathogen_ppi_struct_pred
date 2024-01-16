import os
import pickle
import argparse
import numpy as np
from matplotlib import pyplot as plt


def plot_msa_v2(predictions_dir, sort_lines=False, dpi=300, output_format='png'):
    # adapted from: https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
    with open(os.path.join(predictions_dir, 'features.pkl'), 'rb') as f:
        feature_dict = pickle.load(f)
    seq = feature_dict["msa"][0]
    if "asym_id" in feature_dict:
        Ls = [0]
        k = feature_dict["asym_id"][0]
        for i in feature_dict["asym_id"]:
            if i == k:
                Ls[-1] += 1
            else:
                Ls.append(1)
            k = i
    else:
        Ls = [len(seq)]
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"]

    msa = feature_dict["msa"][:N]
    gap = msa != 21
    qid = msa == seq
    gapid = np.stack([gap[:, Ln[i]:Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:, Ln[i]:Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]
        Nn.append(len(lines_))
        lines.append(lines_)

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)
    plt.figure(figsize=(8, 5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(lines,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
               extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(
        os.path.join(predictions_dir, 'msa_sequence_coverage_plot.{}'.format(output_format)))
    return plt


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Calculate scores, interface statistics and number of seqs in MSAs for AlphaFold2 predicted structures')
    arg_parser.add_argument(
        '-i',
        dest='predictions_dir',
        default=None,
        help='Folder containing the AlphaFold-multimer outputs',
    )
    arg_parser.add_argument(
        '--sort-lines',
        dest='sort_lines',
        action='store_true',
        help='Sorts the sequences by sequence identity to query',
    )
    arg_parser.add_argument(
        '--dpi',
        dest='dpi',
        default=300,
        type=int,
        help='DPI of figure',
    )
    arg_parser.add_argument(
        '--output-format',
        dest='output_format',
        default=None,
        help='Output format of figure',
    )
    args = arg_parser.parse_args()
    plot_msa_v2(**vars(args))
