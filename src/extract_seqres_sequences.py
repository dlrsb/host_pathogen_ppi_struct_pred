import os
import glob
import argparse
from Bio import SeqIO


def pdb_to_fasta(pdb_filepath, output_dir=None):
    records = list(SeqIO.parse(pdb_filepath, "pdb-atom"))
    chains = [records[i].id.split(':')[-1] for i in range(len(records))] # get the chain IDs that are present in this biological assembly file
    print(chains)
    records_to_write = []
    for record in SeqIO.parse(pdb_filepath, "pdb-seqres"): # get seqs from the SEQRES records in the PDB file
        print(record)
        if record.id in chains:
            record.description = record.id
            records_to_write.append(record)
    pdb_dir, filename = os.path.split(pdb_filepath)
    if output_dir is None:
        output_dir = pdb_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    fasta_filepath = os.path.join(output_dir, '{}.fasta'.format(filename.split('.')[0]))
    with open(fasta_filepath, "w") as output_handle:
        SeqIO.write(records_to_write, output_handle, "fasta")


def pdbs_to_fasta(pdb_dir_path, output_dir):
    pdb_paths = glob.glob(os.path.join(pdb_dir_path, '*.pdb'))
    for pdb_path in pdb_paths:
        pdb_to_fasta(pdb_path, output_dir)


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description='Extracts protein sequences from SEQRES records in PDB files and saves them as FASTA files.')
    arg_parser.add_argument(
        '-i',
        dest='pdb_filepath',
        default=None,
        help='Path to the PDB file',
    )
    arg_parser.add_argument(
        '-o',
        dest='output_dir',
        default=None,
        help='Path to output directory',
    )
    args = arg_parser.parse_args()
    pdb_to_fasta(**vars(args))
