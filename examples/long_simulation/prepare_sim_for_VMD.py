import MDAnalysis as mda
from MDAnalysis import transformations as trans
import numpy as np
from MDAnalysis.analysis import align  # Import the align module for inertia tensor calculation

pdb_name = 'system_prep/last_frame_0.pdb'
traj_name = ['traj_50ns_1.dcd', 'traj_50ns_2.dcd', 'traj_50ns_3.dcd']

step = 10

universe_mda = mda.Universe(pdb_name, traj_name)

print('Number of atoms:', universe_mda.atoms.n_atoms)
print('Number of frames:', len(universe_mda.trajectory))
print(f'Using {len(universe_mda.trajectory[::step])} frames')

non_water = universe_mda.select_atoms("not resname HOH WAT SOL TIP3 CL NA")
protein = universe_mda.select_atoms("protein")
system = universe_mda.select_atoms("protein or resname UNK")

transformations = [
    trans.center_in_box(protein),  # center protein
    trans.unwrap(universe_mda.atoms),  # unwrap molecules across PBC (optional, can be removed)
    trans.wrap(universe_mda.atoms),    # wrap back into box (optional, can be removed)
]

print('Applying transformations')
universe_mda.trajectory.add_transformations(*transformations)

print('Writing PDB')
# Write PDB for the first frame
universe_mda.trajectory[0]  # or choose another frame if you prefer
non_water.write("processed_no_water.pdb")

print('Writing DCD')
# --- Write trajectory ---
with mda.Writer("processed_no_water.dcd", non_water.n_atoms) as W:
    for ts in universe_mda.trajectory[::step]:
        W.write(non_water)

print("Done!")