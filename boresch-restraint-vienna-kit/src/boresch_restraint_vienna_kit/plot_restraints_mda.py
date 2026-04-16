import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals

def obtain_bonds(dist_array, universe, a1, a2, step):
  for i, ts in enumerate(universe.trajectory[::step]):
    dist = calc_bonds(
      coords1=a1,
      coords2=a2, 
      box=universe.dimensions
    )
    dist_array[i] = dist[0]
  return dist_array


def obtain_angles(angles_array, universe, a1, a2, a3,step):
  for i, ts in enumerate(universe.trajectory[::step]):
    angle = calc_angles(
      coords1=a1,
      coords2=a2, 
      coords3=a3,
      box=universe.dimensions
      )
    angles_array[i] = angle[0]
  return angles_array


def obtain_dihedrals(dihedral_array, universe, a1, a2, a3, a4,step):
  for i, ts in enumerate(universe.trajectory[::step]):
    dihedral = calc_dihedrals(
      coords1=a1,
      coords2=a2, 
      coords3=a3,
      coords4=a4,
      box=universe.dimensions
      )
    dihedral_array[i] = dihedral[0]
  return dihedral_array

def center_dihedral(array):
  for i in range(len(array)):
    if array[i] < 30:
      array[i] += 360
  return array


def plot_iter(fig, axs,arr1,arr2,arr3,arr4,arr5,arr6,color,label):
  axs[0,0].plot(arr1, color=color, label=label)
  axs[0,1].plot(arr2, color=color )
  axs[1,0].plot(arr3, color=color )
  axs[1,1].plot(arr4, color=color )
  axs[2,0].plot(arr5, color=color )
  axs[2,1].plot(arr6, color=color )


def post_plot(fig, axs,references):
  #! First plot
  axs[0,0].set_title('Distance Host-Guest ($R_A$)')
  axs[0,0].set_ylabel('$R_A$ / Å')
  if references:
    axs[0,0].axhline(8.858, color='k', linestyle='--', label='Reference') 

  #! Second plot 
  axs[0,1].set_title('Angle on the host ($\\theta_A$)')
  axs[0,1].set_ylabel('$\\theta_A$ / $º$')
  if references:
    axs[0,1].axhline(np.degrees(1.234), color='k', linestyle='--', label='Reference') 
  #!Third plot
  axs[1,0].set_title('Angle on the ligand ($\\theta_B$)')
  axs[1,0].set_ylabel('$\\theta_B$ / $º$')
  if references:
    axs[1,0].axhline(np.degrees(1.762), color='k', linestyle='--', label='Reference') 

  #!Fourth plot
  axs[1,1].set_title('Dihedral on the Host ($\\phi_A$)')
  axs[1,1].set_ylabel('$\\phi_A$ / $º$')
  if references:
    axs[1,1].axhline(np.degrees(-0.332), color='k', linestyle='--', label='Reference') 

  #!Fifth plot
  axs[2,0].set_title('Dihedral middle ($\\phi_B$)')
  axs[2,0].set_ylabel('$\\phi_B$ / $º$')
  if references:
    axs[2,0].axhline(np.degrees(0.413), color='k', linestyle='--', label='Reference') 

  #! Sixth plot
  axs[2,1].set_title('Dihedral on the ligand ($\\phi_C$)')
  axs[2,1].set_ylabel('$\\phi_C$ / $º$')
  if references:
    axs[2,1].axhline(np.degrees(-1.802), color='k', linestyle='--', label='Reference') 
  
  axs[0,0].tick_params(labelbottom=False)
  axs[0,1].tick_params(labelbottom=False)

  handles, labels = axs[0,0].get_legend_handles_labels()

  fig.legend(handles, labels,
             loc='center right',
             title='Replicas')

  plt.tight_layout(rect=[0,0,0.90,1])

  # plt.show()


def plot_restraints(pdb_name:str,
                    traj_name,
                    step:int,
                    anchors:list,
                    several_simulations:bool,
                    total_simulations:int,
                    host_idx_corrector:int,
                    guest_idx_corrector:int,
                    references:bool,
                    figure_name:str):
  fig, axs = plt.subplots(3,2,sharex=True, figsize=(14,10))

  colors = plt.cm.tab10(np.linspace(0,1,10))

  for i in range(total_simulations):
    if several_simulations:
      print(f'rep_{i}')
      os.chdir(f'rep_{i}')
    u = mda.Universe(
      pdb_name,
      traj_name
    )

    print('Number of atoms:',u.atoms.n_atoms)
    print('Number of frames:',len(u.trajectory))
    print(f'Using {len(u.trajectory[::step])} frames')

    #* Remember to add one to what you get from the openfe.log
    #! HOST ATOMS !#
    h0 = u.select_atoms(f'id {anchors[0] +host_idx_corrector}') 
    h1 = u.select_atoms(f'id {anchors[1] +host_idx_corrector}') 
    h2 = u.select_atoms(f'id {anchors[2] +host_idx_corrector}') 

    #! GUEST ATOMS !#
    g0 = u.select_atoms(f'id {anchors[3] +guest_idx_corrector}') 
    g1 = u.select_atoms(f'id {anchors[4] +guest_idx_corrector}') 
    g2 = u.select_atoms(f'id {anchors[5] +guest_idx_corrector}') 

    atoms = [h0, h1, h2, g0, g1, g2]
    names = [a[0].name for a in atoms]
    idx = [anchors[i] + (host_idx_corrector if i < 3 else guest_idx_corrector) for i in range(len(anchors))]
    print('\nCheck that your atoms are correct, if not change the correctors in the functions: ')
    print('indexes may be different by some numbers, trust the atom name. \n')
    print(idx)
    print(names)

    #! Computing Parameters !#
    print('\n - Computing Parameters.')
    n_frames = len(u.trajectory[::step])
    
    print('   - Computing distance r_aA.')
    r_A = np.zeros(n_frames)
    r_A = obtain_bonds(r_A, u, h0, g0, step)

    print('   - Computing angle theta_A.')
    theta_A = np.zeros(n_frames)
    theta_A = obtain_angles(theta_A, u, h1, h0, g0, step)
    theta_A = np.degrees(theta_A)

    print('   - Computing angle theta_B.')
    theta_B = np.zeros(n_frames)
    theta_B = obtain_angles(theta_B, u, h0, g0, g1, step)
    theta_B = np.degrees(theta_B)

    print('   - Computing angle phi_A.')
    phi_A = np.zeros(n_frames)
    phi_A = obtain_dihedrals(phi_A, u, h2, h1, h0, g0, step)
    phi_A = np.unwrap(phi_A)

    phi_A = np.degrees(phi_A)

    # phi_A = center_dihedral(phi_A)

    print('   - Computing angle phi_B.')
    phi_B = np.zeros(n_frames)
    phi_B = obtain_dihedrals(phi_B, u, h1, h0, g0, g1, step)
    phi_B = np.unwrap(phi_B)
    phi_B = np.degrees(phi_B)

    # phi_B = center_dihedral(phi_B)

    print('   - Computing angle phi_B.')
    phi_C = np.zeros(n_frames)
    phi_C = obtain_dihedrals(phi_C, u, h0, g0, g1, g2, step)
    phi_C = np.unwrap(phi_C)
    phi_C = np.degrees(phi_C)

    # phi_C = center_dihedral(phi_C)

    print('\n - Plotting.')
    plot_iter(fig, axs,r_A,theta_A,theta_B,phi_A,phi_B,phi_C,color=colors[i],label=f'rep_{i}')
    if several_simulations:
      os.chdir('..')


  post_plot(fig, axs, references)
  plt.savefig(f"{figure_name}.svg", bbox_inches="tight")