import os
import argparse

from score_prediction import score_structure


def score_all(results_dir, native_structures_dir, fasta_files_dir, output_filepath, rename_chains,
              only_best_model=False, select_residues_before_dockq=False, interface_cutoff=10):
    prediction_dirs = [pred_dir for pred_dir in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, pred_dir))]
    for pred_dir in prediction_dirs:
        if fasta_files_dir is not None:
            fasta_path = os.path.join(fasta_files_dir, f'{pred_dir}.fasta')
        else:
            fasta_path = None
        if os.path.exists(os.path.join(results_dir, pred_dir, f'{pred_dir}.done')):
            pred_dir_path = os.path.join(results_dir, pred_dir, pred_dir)

            if native_structures_dir is not None:
                if os.path.exists(os.path.join(native_structures_dir, f'{pred_dir}_fixedmodel.pdb')):
                    native_path = os.path.join(native_structures_dir, f'{pred_dir}_fixedmodel.pdb')
                    if os.path.exists(os.path.join(results_dir, pred_dir, f'{pred_dir}_fixedmodel')):
                        pred_dir_path = os.path.join(results_dir, pred_dir, f'{pred_dir}_fixedmodel')
                else:
                    native_path = os.path.join(native_structures_dir, f'{pred_dir}.pdb')
            else:
                native_path = None

            try:
                print('trying ' + pred_dir)
                score_structure(predictions_dir=pred_dir_path,
                                native_structure_filepath=native_path,
                                input_fasta_filepath=fasta_path,
                                rename_chains=rename_chains,
                                only_best_model=only_best_model,
                                select_residues_before_dockq=select_residues_before_dockq,
                                interface_cutoff=interface_cutoff,
                                output_filepath=output_filepath)
            except Exception as e:
                print('failed ' + pred_dir)
                print(e)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Scores multimer AlphaFold-multimer structures')
    arg_parser.add_argument('results_dir', help='Folder containing Alphafold-multimer output folders for a given dataset')
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
        help='Path to the directory containing the input fasta files',
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
    score_all(**vars(args))
