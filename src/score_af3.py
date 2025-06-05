import os
import shutil
import gzip
import json
import argparse
from Bio.PDB import MMCIFParser, PDBIO
from scoring_metrics import *
from score_prediction import save_results


def score_structure(pair, predictions_dir, native_structure_filepath=None, input_fasta_filepath=None,
                    select_residues_before_dockq=False, interface_cutoff=8, output_filepath=None):
    results = {'name': pair}
    compressed_model_filepath = os.path.join(predictions_dir, f'{pair.lower()}_model.cif.gz')
    results['ranked0_model_path'] = compressed_model_filepath
    model_filepath = os.path.join('/cluster/scratch/dbaptista/score_tmp/', os.path.split(os.path.splitext(compressed_model_filepath)[0])[-1])
    pdb_filepath = f'{os.path.splitext(model_filepath)[0]}.pdb'
    af3_scores_pth = os.path.join(predictions_dir, f'{pair.lower()}_summary_confidences.json.gz')

    # get AF3 scores from *summary_confidences.json
    with gzip.open(af3_scores_pth, 'rt') as f:
        af3_scores = json.load(f)
    af3_scores_selected = {k: v for k, v in af3_scores.items() if k in ['fraction_disordered', 'has_clash', 'iptm', 'ptm', 'ranking_score']}
    results.update(af3_scores_selected)

    # decompress model .cif.gz file
    with gzip.open(compressed_model_filepath, 'rb') as f1:
        with open(model_filepath, 'wb') as f2:
            shutil.copyfileobj(f1, f2)
    # convert .cif to pdb
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(pair, model_filepath)
    io = PDBIO()
    io.set_structure(struct)
    io.save(pdb_filepath)

    # calculate pDockQ as defined in the "Towards a structurally resolved human protein interaction network" paper,
    # with a user-defined cutoff value
    pdockq_huintaf2_results = pdockq_huintaf2_score(pdb_filepath, cutoff=interface_cutoff)
    results.update(pdockq_huintaf2_results)

    # Calculate MM-align TM-score and DockQ score
    if native_structure_filepath is not None:
        # MM-align TM-score
        tmscore_results = tm_score(native_structure_filepath=native_structure_filepath,
                                   predicted_structure_filepath=pdb_filepath)
        results.update(tmscore_results)
        # DockQ
        dockq_results = dockq_score(native_structure_filepath=native_structure_filepath,
                                    predicted_structure_filepath=pdb_filepath,
                                    fasta_filepath=input_fasta_filepath,
                                    select_residues=select_residues_before_dockq)
        results.update(dockq_results)

    if output_filepath is None:
        output_filepath = os.path.join(predictions_dir, 'results.csv')
    save_results(results_dict=results, output_filepath=output_filepath)

    os.remove(model_filepath)
    os.remove(pdb_filepath)


def score_all(results_dir, native_structures_dir, fasta_files_dir, output_filepath,
              select_residues_before_dockq=False, interface_cutoff=8):
    # fasta_files_dir - directory with the AF3 input sequences, which will be the fullseq_fasta_files used for AF2
    prediction_dirs = [pred_dir for pred_dir in os.listdir(results_dir) if
                       os.path.isdir(os.path.join(results_dir, pred_dir))]

    for pred_dir in prediction_dirs:
        pair_name = pred_dir.upper()
        if fasta_files_dir is not None:
            fasta_filepath = os.path.join(fasta_files_dir, f'{pair_name}.fasta')
        else:
            fasta_filepath = None

        if native_structures_dir is not None:
            if os.path.exists(os.path.join(native_structures_dir, f'{pair_name}_fixedmodel.pdb')):
                native_filepath = os.path.join(native_structures_dir, f'{pair_name}_fixedmodel.pdb')
            else:
                native_filepath = os.path.join(native_structures_dir, f'{pair_name}.pdb')
        else:
            native_filepath = None

        try:
            score_structure(pair_name, os.path.join(results_dir, pred_dir), native_filepath, fasta_filepath,
                            select_residues_before_dockq=select_residues_before_dockq,
                            interface_cutoff=interface_cutoff, output_filepath=output_filepath)
        except Exception as e:
            print('failed ' + pair_name)
            print(e)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Scores AF3 structures')
    arg_parser.add_argument('results_dir', help='Folder containing Alphafold3 output folders for a given dataset')
    arg_parser.add_argument(
        '-n',
        dest='native_structures_dir',
        default=None,
        help='Path to the directory containing the native (experimentally-determined) PDB files',
    )
    arg_parser.add_argument(
        '-f',
        dest='fasta_files_dir',
        default=None,
        help='Path to the directory containing the input fasta files (the ones used for AF2 and that were used to create the input JSON files for AF3)',
    )
    arg_parser.add_argument(
        '-o',
        dest='output_filepath',
        default=None,
        help='Path to output file (results CSV)',
    )
    arg_parser.add_argument(
        '--interface-cutoff',
        dest='interface_cutoff',
        default=8,
        type=int,
        help='Distance cutoff value to define the interface',
    )
    arg_parser.add_argument(
        '--select-residues-before-dockq',
        dest='select_residues_before_dockq',
        action='store_true',
        help='Only calculate DockQ based on the residues that also exist in the native file',
    )
    args = arg_parser.parse_args()
    score_all(**vars(args))
