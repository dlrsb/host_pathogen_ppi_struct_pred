import os
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


def rename_protein_chains(pdb_filepath, fasta_filepath):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    name = os.path.split(pdb_filepath)[-1].split('.')[0]
    structure = parser.get_structure(name, pdb_filepath)

    fasta_records = list(SeqIO.parse(fasta_filepath, 'fasta'))

    blosum62_alphabet = substitution_matrices.load('BLOSUM62').alphabet
    for record in fasta_records:
        seq = record.seq
        for char in seq:
            if char not in blosum62_alphabet:
                seq = seq.replace(char, 'X') # AF2 will also use 'X' instead of chars like 'U'
        record.seq = seq

    model_records = list(SeqIO.parse(pdb_filepath, 'pdb-atom'))
    chain_id_map = {}
    for record in model_records:
        aln_scores = {}
        for fasta_record in fasta_records:
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
            aligner.open_gap_score = -10.0
            aligner.extend_gap_score = -0.5
            aln_scores[fasta_record.id] = aligner.score(record.seq, fasta_record.seq)
        print(aln_scores)
        chain_id_map[record.id.split(':')[-1]] = max(aln_scores, key=lambda key: aln_scores[key])
    chain_id_map_temp_ids = [{k: v.lower() for k, v in chain_id_map.items()}, {v.lower(): v for k, v in chain_id_map.items()}]
    print(chain_id_map)
    print(chain_id_map_temp_ids)
    for map in chain_id_map_temp_ids:
        for model in structure:
            for chain in model:
                if chain.id in map:
                    chain.id = map[chain.id]

    io.set_structure(structure)
    pdb_dir, filename = os.path.split(pdb_filepath)
    new_filename = '{}_chains_renamed.pdb'.format(filename.split('.')[0])
    io.save(os.path.join(pdb_dir, new_filename))
