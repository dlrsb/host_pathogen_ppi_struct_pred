import os
import argparse
from Bio.PDB import PDBIO
from utils import load_structure


def rewrite_pdb(pdb_filepath, output_filepath):
    struct = load_structure(pdb_filepath)
    io = PDBIO()
    io.set_structure(struct)
    io.save(output_filepath)
    return output_filepath


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Rewrites PDB file')
    arg_parser.add_argument(
        '-f',
        dest='pdb_filepath',
        default=None,
        help='Path to the PDB file to be rewritten',
    )
    arg_parser.add_argument(
        '-o',
        dest='output_filepath',
        default=None,
        help='Path to the output file',
    )
    args = arg_parser.parse_args()
    rewrite_pdb(**vars(args))