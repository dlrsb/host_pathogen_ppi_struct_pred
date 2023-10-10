#!/usr/bin/env python
"""
Auxiliary functions to compute the pDockQ score for the 3rd week's practical of
the course 551-1299-00L Bioinformatics HS2022 @ ETHZ

Uses BioPython objects.

__author__ = "Miguel Correa Marrero"
"""

import numpy as np
from warnings import warn


def get_mean_interface_plddt(structure, interface_a, interface_b):
    """
    Calculate the mean pLDDT of residues found at the interface.

    Arguments
    ---------
    structure: Bio.PDB.Structure.Structure object
    interface_a: array-like, list of residue indexes at the interface in chain A
    interface_b: array-like, list of residue indexes at the interface in chain B

    Returns
    -------
    mean_interface_plddt: float, mean pLDDT at interface

    """
    if len(interface_a) > 0 and len(interface_b) > 0:
        plddts_a = get_chain_plddts(structure["A"], interface_a)
        plddts_b = get_chain_plddts(structure["B"], interface_b)
        mean_interface_plddt = np.mean(np.concatenate([plddts_a, plddts_b]))
    else:
        warn("No interface residues; returning NaN")
        mean_interface_plddt = np.nan

    return mean_interface_plddt


def get_chain_plddts(chain, residue_idxs):
    """
    Retrieve the pLDDT for each specified residue in a given protein chain.
    The pLDDT is stored in the B-factor column in the PDB file

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object
    residue_idxs: array-like, list of residue indexes to retrieve pLDDT for

    Returns
    -------
    chain_plddts: array-like, pLDDTs for each residue
    """
    chain_plddts = np.zeros(len(residue_idxs))
    for idx, residue_idx in np.ndenumerate(residue_idxs):
        chain_plddts[idx] = chain[int(residue_idx)]["CA"].get_bfactor()
    return chain_plddts


def get_number_interface_residues(interface_a, interface_b):
    """
    Calculates the total number of residues found at the interface
    """
    return len(interface_a) + len(interface_b)


def calculate_pdockq(avg_interface_plddt, no_interface_residues):
    x = avg_interface_plddt * np.log(no_interface_residues)
    exp_term = np.exp(-0.03148 * (x - 388.06))
    denom = 1 + exp_term
    pdockq = (0.707 / denom) + 0.03138
    assert not (pdockq < 0)
    assert not (pdockq > 1)
    return pdockq
