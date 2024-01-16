import os
from Bio import SeqIO


def create_single_prot_fasta_files(pairs_fasta_dir_path, individual_fasta_dir_path):
    fasta_files = os.listdir(pairs_fasta_dir_path)
    for fasta_file in fasta_files:
        for record in SeqIO.parse(os.path.join(pairs_fasta_dir_path, fasta_file), "fasta"):
            print(record.id)
            output_fasta_path = os.path.join(individual_fasta_dir_path, '{}.fasta'.format(record.id))
            if not os.path.exists(output_fasta_path):
                SeqIO.write(record, output_fasta_path, 'fasta-2line')


if __name__ == '__main__':
    create_single_prot_fasta_files('../data/sars_mers_sars2_human_ppi/protein_pairs',
                                   '../data/sars_mers_sars2_human_ppi/individual_proteins')