#!/usr/bin/env python3

import sys
import numpy as np
import argparse
import re
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
#import yaml


# Ny version
# m = np.sqrt(np.sum((a[:,np.newaxis,:] -a[np.newaxis,:,:])**2, axis=2))
# n=np.where(m<=6.0)

# mat = np.append(chi_coords,chj_coords,axis=0)
# a_min_b = (mat[:,np.newaxis,:] -mat[np.newaxis,:,:])
# dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0))
# contact_dists = dists[l1:,:l1]
# t=8
# contacts = np.argwhere(contact_dists<=t)


def pDockQ(structure, cutoff, disocut1=50, disocut2=70, disocut3=90):
    i = 0
    tiny = 1.e-20
    chains = []
    for chain in structure.get_chains():
        chains += [chain]
        i += 1
        # print (chains)
    interface_residues1 = []
    interface_residues2 = []
    for res1 in chains[0]:
        for res2 in chains[1]:
            # Atom-atom distance
            # print (res1,res2)
            test = False
            for i in res1:
                if test: break
                for j in res2:
                    dist = np.linalg.norm(i.coord - j.coord)
                    # scipy.spatial.distance.euclidian
                    # dist = distance.euclidean(coords[i], coords[len(ch1_res_nos)+j]) #Need to add l1 to get the right coords
                    if dist < cutoff:
                        # Save residues
                        hetflag, resseq, icode = res1.get_id()
                        interface_residues1.append(resseq)
                        hetflag, resseq, icode = res2.get_id()
                        interface_residues2.append(resseq)
                        test = True
                        break
                    elif dist > 2 * cutoff:  # To speed up things
                        test = True
                        break

    # print (interface_residues1,interface_residues2)
    # yaml.dump({
    #     'resseq1': [int(a) for a in np.unique(interface_residues1)],
    #     'resseq2': [int(a) for a in np.unique(interface_residues2)],
    # }, stream=sys.stdout, default_flow_style=None)
    resseq1 = [int(a) for a in np.unique(interface_residues1)]
    resseq2 = [int(a) for a in np.unique(interface_residues2)]
    # interface_res_num.append(np.unique(interface_residues).shape[0])
    # atoms, residue_numbers, coords = np.array(atoms), np.array(residue_numbers), np.array(coords)
    if1 = np.unique(interface_residues1)
    if2 = np.unique(interface_residues2)
    NumRes = if1.shape[0] + if2.shape[0]
    i = tiny
    b = 0
    b1 = 0
    b2 = 0
    i1 = 0
    i2 = 0
    NumDiso1 = [0, 0, 0, 0]
    NumDiso2 = [0, 0, 0, 0]
    for res in chains[0]:
        b1 += res['CA'].get_bfactor()
        i1 += 1
        if res['CA'].get_bfactor() > disocut3:  # >90
            NumDiso1[0] += 1
        elif res['CA'].get_bfactor() > disocut2:  # 70-90
            NumDiso1[1] += 1
        elif res['CA'].get_bfactor() > disocut1:  # 50-70
            NumDiso1[2] += 1
        else:  # <50
            NumDiso1[3] += 1
        if res.id[1] in if1:
            b += res['CA'].get_bfactor()
            i += 1
    for res in chains[1]:
        b2 += res['CA'].get_bfactor()
        i2 += 1
        if res['CA'].get_bfactor() > disocut3:  # >90
            NumDiso2[0] += 1
        elif res['CA'].get_bfactor() > disocut2:  # 70-90
            NumDiso2[1] += 1
        elif res['CA'].get_bfactor() > disocut1:  # 50-70
            NumDiso2[2] += 1
        else:  # <50
            NumDiso2[3] += 1
        if res.id[1] in if2:
            b += res['CA'].get_bfactor()
            i += 1
    IF_plDDT = b / i
    plDDT1 = b1 / i1
    plDDT2 = b2 / i2
    # print ("test",b,b1,b2,i,i1,i2,NumDiso1,NumDiso2)
    # Get res nos
    # Get chain cut
    # ch1_res_nos = np.argwhere(residue_numbers<=l1)[:,0] #All residue numbers
    # ch2_res_nos =  np.argwhere(residue_numbers>l1)[:,0]

    # print (NumRes,IF_plDDT)
    return (NumRes, IF_plDDT, plDDT1, plDDT2, NumDiso1, NumDiso2, i1, i2, resseq1, resseq2)


def sigmoid(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return (y)


# arg_parser = argparse.ArgumentParser(description="Calculates pDockQ from NumRes and IF_plDDT")
# group = arg_parser.add_mutually_exclusive_group(required=True)
# group.add_argument("-p", "--pdb", type=argparse.FileType('r'), help="Input pdb file pLddt values in bfactor columns")
# group.add_argument("-c", "--cif", type=argparse.FileType('r'), help="Input cif file plddt values in bfactor columns")
# group.add_argument("-f", "--file", type=argparse.FileType('r'), help="Input precalculated pLDDT file.")
# group.add_argument("-n", "--NumRes", type=float, required=False, help="Number of residues in interface")
#
# arg_parser.add_argument("-F", "--factor", type=float, required=False, default=1.0,
#                         help="Factor to multiply B-factor for (sometimes 100)")
# arg_parser.add_argument("-C", "--cutoff", type=float, required=False, default=10.0,
#                         help="Cutoff for defining distances")
#
# arg_parser.add_argument("-i", "--IF_plDDT", type=float, required=False, help="Average plDDT in interface residues")
# arg_parser.add_argument("-v", "--verbose", action='store_true', required=False,
#                         help="Average plDDT in interface residues")
#
# args = arg_parser.parse_args()
# cutoff = args.cutoff
#
# # popt=[7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]
# # popt2=[-6.13932499e-01,  3.87553637e+02, -1.82925395e-02,  6.32058555e-01] # optimized on DockQall from human
# # Testing using different cutoffs.
# if cutoff <= 5:
#     # 5 SpearmanrResult(correlation=0.7647585237390458, pvalue=4.749030057232305e-280)
#     popt = [6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
#     # 0.7805034405869632
# elif cutoff <= 6:
#     # 6 SpearmanrResult(correlation=0.7708834427476546, pvalue=2.707297682746201e-287)
#     popt = [7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
#     # 0.7871982094514278
# elif cutoff <= 7:
#     # 7 SpearmanrResult(correlation=0.7709518988131879, pvalue=2.2402500804327052e-287)
#     popt = [7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
#     # 0.7859609807320201
# elif cutoff <= 8:
#     # 8 SpearmanrResult(correlation=0.7632969367380509, pvalue=2.3583905451705336e-278)
#     popt = [7.18442739e-01, 3.60791204e+02, 3.01635944e-02, 2.04076969e-02]
#     # 0.7764648775754815
# elif cutoff <= 9:
#     # 9 SpearmanrResult(correlation=0.7496303495195178, pvalue=4.539049646719674e-263)
#     popt = [7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
#     # 0.7608417399783565
# elif cutoff <= 10:
#     # 10 SpearmanrResult(correlation=0.7330653937901442, pvalue=7.988440779428826e-246)
#     # popt=[7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02] # Newly optimizes
#     popt = [7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]  # used in Interaction studies.
#     # 0.7431426093979494
# elif cutoff <= 11:
#     # 11 SpearmanrResult(correlation=0.71288058226417, pvalue=1.7542846392453894e-226)
#     popt = [7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
#     # 0.7219615906164123
# else:
#     # 12 SpearmanrResult(correlation=0.6938911161134763, pvalue=9.284495013784153e-210)
#     popt = [7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]
#     # 0.7023000652310362
#
# cutoff2 = 3
# Numoverlap = 0
# overlap = 0
# Name = ""
# if (args.file):
#     # file = open(args.file, 'r')
#     file = args.file
#     lines = file.readlines()
#     for line in lines:
#         temp = line.split()
#         if (temp[0] == "IF_NumRes:"):
#             NumRes = float(temp[1])
#         if (temp[0] == "IF_pLDDT"):
#             IF_plDDT = float(temp[1])
#     Name = args.file.name
# elif (args.cif):
#     name = args.cif.name
#     bio_parser = MMCIFParser()
#     structure_file = args.cif
#     structure_id = args.cif.name[:-4]
#     structure = bio_parser.get_structure(structure_id, structure_file)
#     NumRes, IF_plDDT, plDDT1, plDDT2, Diso1, Diso2, len1, len2, residues_seq1, residues_seq2 = pDockQ(structure, cutoff)
#     if (args.verbose):
#         NumResOverlap, IF_plDDTOverlap, plDDT1Overlap, plDDT2overlap, Diso1overlap, Diso2overlap, len1, len2, \
#             residues_seq1_overlap, residues_seq2_overlap = pDockQ(structure, cutoff2)
#
#
# elif (args.pdb):
#     Name = args.pdb.name
#     bio_parser = PDBParser()
#     structure_file = args.pdb
#     structure_id = args.pdb.name[:-4]
#     structure = bio_parser.get_structure(structure_id, structure_file)
#     NumRes, IF_plDDT, plDDT1, plDDT2, Diso1, Diso2, len1, len2, residues_seq1, residues_seq2= pDockQ(structure, cutoff)
#     if (args.verbose):
#         NumResOverlap, IF_plDDTOverlap, plDDT1Overlap, plDDT2overlap, Diso1overlap, Diso2overlap, len1, len2,
#         residues_seq1_overlap, residues_seq2_overlap = pDockQ(
#             structure, cutoff2)
# else:
#     NumRes = args.NumRes
#     IF_plDDT = args.IF_plDDT
#
# tiny = 1.e-20
# name = re.sub(r'.*/', '', Name)
# name = re.sub(r'.pdb$', '', name)
# name = re.sub(r'.cif$', '', name)
# dir = re.sub(r'/[\w_0-9A-Za-z]+\.\w+', '', Name)
# dir = re.sub(r'.*/', '', dir)
#
# # print (NumRes,tiny,IF_plDDT, popt)
# pDockQ = sigmoid(np.log(NumRes + tiny) * IF_plDDT * args.factor, *popt)
#
# if (args.verbose):
#     # print ("PdockQ: %.3f\nNumRes: %d\nIF_plDDT: %.2f\n" % (pDockQ,NumRes,IF_plDDT))
#     print(
#         "Name,Dir,pDockQ,NumRes,IF_plDDT,plDDT1,plDDT2,NumDiso1+90,NumDiso1-70-90,NumDiso1-50-70,NumDiso1-50,NumDiso2+90,NumDiso2-70-90,NumDiso2-50-70,NumDiso2-50,NumOverlap,len1,len2")
#     print("%s,%s,%f,%d,%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d" % (
#     name, dir, pDockQ, NumRes, IF_plDDT, plDDT1, plDDT2, Diso1[0], Diso1[1], Diso1[2], Diso1[3], Diso2[0], Diso2[1],
#     Diso2[2], Diso2[3], NumResOverlap, len1, len2))
# else:
#     # print (np.round(pDockQ,3))
#     yaml.dump({
#         'pDockQ': round(float(pDockQ), 3),
#     }, stream=sys.stdout)
