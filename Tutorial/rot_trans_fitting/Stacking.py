#!/usr/bin/env python
# coding: utf-8

# # Rotation and translation fitting of Amyloid-B

import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.core.universe import Merge
import numpy as np
from copy import deepcopy
import string


raw_structure_path = '2beg.pdb'
universe = mda.Universe(raw_structure_path)

center_universe = universe.segments[2:4]
chainA = universe.segments[2].atoms
chainB = universe.segments[3].atoms

# Write the first chain we use in the stacking, so we can 
chainA.write('chainA.pdb')

# Only using the C-alphas for fitting
ref_atoms = chainB.select_atoms('protein and name CA').positions
mob_atoms = chainA.select_atoms('protein and name CA').positions

# Get the center of mass coordinates
ref_COG = np.mean(ref_atoms, axis=0)
mob_COG = np.mean(mob_atoms, axis=0)

# Get coordinates and rebase the center of mass to 0,0,0
ref_origin = ref_atoms - ref_COG
mob_origin = mob_atoms - mob_COG

# Compute the rotation matrix and RMSD (rotation around COM leaves COM unchanged)
R, rmsd = rotation_matrix(mob_origin, ref_origin)

# Returns the fitted positions for only the CA
mob_aligned = (mob_atoms - mob_COG) @ R + ref_COG

# Better yet is to generate a new empty universe in memory without passing by the HD
new_chain = mda.Universe('chainA.pdb').atoms 
# Apply the rot+trans operation to the new copy
new_chain.positions = (chainA.positions - mob_COG) @ R + ref_COG # All atoms from chainA projected on chainB


# Combine the original A with the new_chain to use as a SnakeOil template

# THIS IS NASTY AND WE LIKE IT THAT WAY!
# Set segment names
#chainA.segments.segids = ['A'] 
#new_chain.segments.segids = ['B']

# Set the chainIDs correctly, martinize depends on this!
# chainA.chainIDs = ['A'] * len(chainA.atoms)
# new_chain.chainIDs = ['B'] * len(chainA.atoms)
# combined_universe = Merge(chainA, new_chain)
# combined_universe.atoms.write('stacked_chains.pdb')


# Make a fiber with 26 (all letters of the alphabet as chains) monomers stacked ontop of eachother, where the stacking is based on the R (rotation matrix)
# and based on the translation matricces (ref_COG, mob_COG)

alphabet = list(string.ascii_uppercase)

initial_chain = mda.Universe('chainA.pdb').atoms 
chainA.chainIDs = ['A'] * len(chainA.atoms)
old_chain = None

for index, chainid in enumerate(alphabet):
    if index == 0:
        old_chain = initial_chain
        old_chain.chainIDs = [chainid] * len(chainA.atoms)
    elif index < 26:
        new_chain = mda.Universe('chainA.pdb').atoms 
        new_chain.positions = (old_chain.positions - mob_COG) @ R + ref_COG
        new_chain.chainIDs = [chainid] * len(chainA.atoms)
        if index == 1:
            combined_universe = Merge(chainA, new_chain)
            combined_universe.atoms.write('multiple_stacked_chains.pdb')
            old_chain = new_chain
        else:
            combined_universe = Merge(combined_universe.atoms, new_chain)
            combined_universe.atoms.write('multiple_stacked_chains.pdb')
            old_chain = new_chain

