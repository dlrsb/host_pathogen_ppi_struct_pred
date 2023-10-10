import os
import json
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


# def rename_protein_chains(pdb_filepath):
#     prediction_dir = os.path.split(pdb_filepath)[0]
#     with open(
#             os.path.join(prediction_dir, 'msas', 'chain_id_map.json')) as f:
#         chain_id_json = json.load(f)
#     chain_id_map = {k: chain_id_json[k]['description'] for k in chain_id_json}
#     chain_id_map_temp_ids = [{k: v.lower() for k, v in chain_id_map.items()}, {v.lower(): v for k, v in chain_id_map.items()}]
#     #print(chain_id_map_temp_ids)
#
#     parser = PDBParser(QUIET=True)
#     io = PDBIO()
#     name = os.path.split(pdb_filepath)[-1].split('.')[0]
#     structure = parser.get_structure(name, pdb_filepath)
#     # for model in structure:
#     #     for chain in model:
#     #         chain.id = chain_id_map[chain.id]
#     try:
#         for map in chain_id_map_temp_ids: # need to do this when I want to change a chain ID to an ID that already exists in the original file
#             for model in structure:
#                 for chain in model:
#                     if chain.id in map:
#                         chain.id = map[chain.id]
#     except ValueError as err: # TODO: test if this approach always works well and if it does, see if I can just use this part for all cases instead
#         #print(err)
#         model_records = list(SeqIO.parse(pdb_filepath, "pdb-atom"))
#         #print(len(model_records))
#         chain_id_map2 = {}
#         for record in model_records:
#             aln_scores = {}
#             for id in chain_id_json:
#                 # if record.seq == fasta_seq:
#                 #     chain_id_map2[record.id.split(':')[-1]] = original_id
#                 fasta_seq = chain_id_json[id]['sequence']
#                 original_id = chain_id_json[id]['description']
#                 aligner = PairwiseAligner()
#                 aligner.mode = 'global'
#                 aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
#                 aligner.open_gap_score = -10.0
#                 aligner.extend_gap_score = -0.5
#                 aln_scores[original_id] = aligner.score(record.seq, fasta_seq)
#             chain_id_map2[record.id.split(':')[-1]] = max(aln_scores, key=lambda key: aln_scores[key])
#         chain_id_map_temp_ids2 = [{k: v.lower() for k, v in chain_id_map2.items()}, {v.lower(): v for k, v in chain_id_map2.items()}]
#         print(chain_id_map2)
#         print(chain_id_map_temp_ids2)
#         structure = parser.get_structure(name, pdb_filepath)
#         for map in chain_id_map_temp_ids2:
#             for model in structure:
#                 for chain in model:
#                     if chain.id in map:
#                         chain.id = map[chain.id]
#
#     io.set_structure(structure)
#     pdb_dir, filename = os.path.split(pdb_filepath)
#     new_filename = '{}_chains_renamed.pdb'.format(filename.split('.')[0])
#     io.save(os.path.join(pdb_dir, new_filename))


def rename_protein_chains(pdb_filepath):
    prediction_dir = os.path.split(pdb_filepath)[0]
    with open(
            os.path.join(prediction_dir, 'msas', 'chain_id_map.json')) as f:
        chain_id_json = json.load(f)
    print(chain_id_json)

    parser = PDBParser(QUIET=True)
    io = PDBIO()
    name = os.path.split(pdb_filepath)[-1].split('.')[0]
    structure = parser.get_structure(name, pdb_filepath)

    model_records = list(SeqIO.parse(pdb_filepath, "pdb-atom"))
    chain_id_map = {}
    for record in model_records:
        aln_scores = {}
        for id in chain_id_json:
            fasta_seq = chain_id_json[id]['sequence']
            original_id = chain_id_json[id]['description']
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
            aligner.open_gap_score = -10.0
            aligner.extend_gap_score = -0.5
            aln_scores[original_id] = aligner.score(record.seq, fasta_seq)
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