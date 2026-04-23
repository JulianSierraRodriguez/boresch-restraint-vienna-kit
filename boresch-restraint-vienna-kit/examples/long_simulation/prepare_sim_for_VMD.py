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

# print('Computing inertia tensor')
# # --- Compute the moment of inertia tensor of the protein ---
# inertia_tensor = protein.moment_of_inertia(wrap=False, unwrap=False, compound='group')
# eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
# rotation_matrix = eigenvectors.T

# # Compute the principal axes (Eigenvectors of the inertia tensor)
# print(eigenvalues)
# quit()
# # The eigenvectors are the principal axes, which we can use to align the protein
#   # Transpose to get the rotation matrix

# # Remove transformations like unwrap and wrap if not necessary for alignment
transformations = [
    trans.center_in_box(protein),  # center protein
    trans.unwrap(universe_mda.atoms),  # unwrap molecules across PBC (optional, can be removed)
    trans.wrap(universe_mda.atoms),    # wrap back into box (optional, can be removed)
]

print('Applying transformations')
universe_mda.trajectory.add_transformations(*transformations)

# def apply_rotation(universe, step):
#     """Rotate the atom coordinates using the given rotation matrix and return the updated universe."""
#     for ts in universe.trajectory[::step]:
#         inertia_tensor = system.moment_of_inertia(wrap=False, unwrap=False, compound='group')
#         eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
#         rotation_matrix = eigenvectors.T
#         # Get the current coordinates
#         coords = system.positions  # The current protein atom coordinates
#         # Apply the rotation matrix: new_coords = coords @ rotation_matrix
#         rotated_coords = coords.dot(rotation_matrix)
#         # Update the atom positions with the rotated coordinates
#         system.positions = rotated_coords
#     return universe  # Returning the updated universe with rotated coordinates

# print('Applying rotation')
# # Call the function and return the updated universe
# updated_universe = apply_rotation(universe_mda, step)

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