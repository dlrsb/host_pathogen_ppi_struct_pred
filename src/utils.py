import os
import glob
import pickle
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser


def load_structure(pdb_filepath):
    parser = PDBParser(QUIET=True)
    name = os.path.split(pdb_filepath)[-1].split('.')[0]
    structure = parser.get_structure(name, pdb_filepath)[0] # gets first model
    return structure


def align_seqs(structure1_filepath, structure2_filepath):
    # based on https://gist.github.com/JoaoRodrigues/e3a4f2139d10888c679eb1657a4d7080
    # get sequences
    struct1_seqs = SeqIO.to_dict(SeqIO.parse(structure1_filepath, 'pdb-atom'))
    struct2_seqs = SeqIO.to_dict(SeqIO.parse(structure2_filepath, 'pdb-atom'))

    # align chains based on id (need to rename chains in 1 of the files beforehand)
    chains = struct1_seqs.keys()
    alignments = {}
    for chain in chains:
        seq1 = struct1_seqs[chain].seq.replace('X', '')
        seq2 = struct2_seqs[chain].seq.replace('X', '')
        # TODO: Switch to Bio.Align.PairwiseAligner as pairwise2 is deprecated.
        #  https://biopython.org/docs/1.76/api/Bio.Align.html
        alignment = pairwise2.align.globalds(
            seq1,
            seq2,
            substitution_matrices.load('BLOSUM62'),
            one_alignment_only=True,
            open=-10.0,
            extend=-0.5,
            penalize_end_gaps=(False, False),
        )
        alignments[chain.split(':')[-1]] = alignment[0]

    return alignments


def reduce_results_file(results_pkl_filepath):
    scores_to_keep = ['num_recycles', 'predicted_aligned_error', 'plddt', 'ptm', 'iptm', 'ranking_confidence']
    # scores_to_keep = ['num_recycles', 'ptm', 'iptm', 'ranking_confidence']
    with open(results_pkl_filepath, 'rb') as f:
        af_results = pickle.load(f)
    selected_scores = {k: v for k, v in af_results.items() if k in scores_to_keep}
    new_pth = '{}_reduced{}'.format(os.path.splitext(results_pkl_filepath)[0], os.path.splitext(results_pkl_filepath)[1])
    with open(new_pth, 'wb') as f:
        pickle.dump(selected_scores, f)
    os.remove(results_pkl_filepath)
