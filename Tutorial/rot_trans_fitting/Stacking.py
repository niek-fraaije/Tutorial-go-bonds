#!/usr/bin/env python
# coding: utf-8

# # Rotation and translation fitting of Amyloid-B

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.core.universe import Merge
import numpy as np

raw_structure_path = '2beg.pdb'
universe = mda.Universe(raw_structure_path)

# Extract atoms
atomsA = universe.segments[2].atoms
atomsB = universe.segments[3].atoms

# Get positions and center
posA = atomsA.positions - atomsA.center_of_mass()
posB = atomsB.positions - atomsB.center_of_mass()

# Compute rotation matrix R to align A to B
R, rmsd = align.rotation_matrix(posA, posB)
T = atomsB.center_of_mass() - atomsA.center_of_mass()  # translation A → B

# Prepare list of positions
positions = []
n_copies = 26

# Initialize with chainA positions
current_pos = atomsA.positions.copy()
current_com = atomsA.center_of_mass()
positions.append(current_pos)

for _ in range(1, n_copies):
    # Apply rotation and translation
    rotated = np.dot(current_pos - current_com, R.T)
    new_com = current_com + T
    new_pos = rotated + new_com

    positions.append(new_pos)
    current_pos = new_pos
    current_com = new_com

# --- Build a new Universe with all 26 chains ---

# Create one merged topology from chainA * 26

stacked = Merge(*[atomsA.copy() for _ in range(n_copies)])

# Set positions
for i, coords in enumerate(positions):
    stacked.segments[i].atoms.positions = coords
    # Set chainID (e.g., A, B, ..., Z, AA, AB, ...)
    if i < 26:
        chain_id = chr(65 + i)  # A–Z
    else:
        chain_id = chr(65 + (i // 26) - 1) + chr(65 + (i % 26))  # AA, AB, ...
    stacked.segments[i].atoms.chainIDs = chain_id


# --- Write output ---
stacked.atoms.write("multiple_stacked_chains.pdb")
