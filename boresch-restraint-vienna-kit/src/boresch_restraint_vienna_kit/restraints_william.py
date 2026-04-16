'''
Procedure obtained from Zhiyi Wu et al. (J. Chem. Theory Comput. (2025)) and coded by Julian Sierra Rodriguez.
'''

#? Imports
import numpy as np
from rdkit import Chem
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals
from rich.pretty import Pretty
from rich.console import Console
from scipy.stats import circvar

#? constants (optional)
global console

console = Console(width=200)

#? core functions

def find_guest_candidates(guest_sdf_name,debug_info):
  #! Step 1: Finding Guest candidates !#

  print('\n 1. Looking for guest anchor candidates.\n')

  # Obtaining the guest molecule in rdkit.
  supp = Chem.SDMolSupplier(guest_sdf_name, removeHs=False)
  guest_rdmol = supp[0]

  # We search the guest molecule for atoms connected to at least two heavy atoms, for example the C in a methyl group would not be included as a candidate.
  
  guest_candidates = [atom for atom in guest_rdmol.GetAtoms() if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]) >= 2]


  # From the array with the guest candidates we are generating an array with the atom indexes from the candidates. 
  guest_candidates_ids  = [atom.GetIdx()+2 for atom in guest_candidates]

  if debug_info:
    # Additional information: We compute the number of atoms of each type selected as candidates.
    c,n,o,other = 0,0,0,0
    other_symb = []
    for atom in guest_candidates:
      symb = atom.GetSymbol()
      if symb == 'C':
        c += 1
      elif symb == 'N':
        n += 1
      elif symb == 'O':
        o += 1
      else:
        other += 1
        other_symb.append(symb)

    print('AMOUNTS OF ATOM CANDIDATES IN GUEST')
    print(f'C = {c}')
    print(f'N = {n}')
    print(f'O = {o}')
    print(f'Other = {other}')
    if other_symb != []:
      print('Other atoms:',other_symb)
    print(f'total = {c+n+o+other}\n')

  print('Indexes of the guest\'s candidate atoms (only guest):')
  print(guest_candidates_ids,'\n')
  return guest_candidates_ids

def setting_up_mda(pdb_name,traj_name,guest_resname,guest_candidates_ids):
  #! MDAnalysis set up !#

  #* Setting up the MDAnalysis universe*#

  u = mda.Universe(
    pdb_name,
    traj_name
  )

  u.atoms.guess_bonds(
    vdwradii={'Na' : 0,
              'Cl' : 0,
              'Br': 1.85}
  )

  #* Coverting guest candidate list to MDA indexes*#

  mda_guest_candidates = []
  # mda_guest_candidates_elem = []
  mda_guest_candidates_idx = []
  ligand_atoms = u.select_atoms(f'resname {guest_resname}')
  for i in guest_candidates_ids:
    atom = ligand_atoms[i]
    mda_guest_candidates.append(atom)
    # mda_guest_candidates_elem.append(atom.element)
    mda_guest_candidates_idx.append(int(atom.index))
  
  

  print('\nIndexes of the guest\'s candidate atoms (full system):')
  print(mda_guest_candidates_idx,'\n')
  mda_guest_candidates_names = [
      u.atoms[i-2].name
      for i in mda_guest_candidates_idx
  ]  
  print(mda_guest_candidates_names)
  return mda_guest_candidates_idx, u, ligand_atoms

def hbond_filter1(universe_mda,ligand_atoms, hbonds, debug_info):
  #* H-bond Filter 1: Generating array where the donnor or acceptor atom is in the guest. *#

  print(f'\nFilter 1: donnor or acceptor atom must be in the guest.')

  protein_atoms = universe_mda.select_atoms("protein")

  host_guest_hbonds = []

  if debug_info:
    print('\n Human readable array from HydrogenBondAnalysis, MDAnalysis:\n')

  for h in hbonds:
    frame, donor, hydrogen, acceptor, dist, angle = h
    donor_atom = universe_mda.atoms[int(donor)]
    acceptor_atom = universe_mda.atoms[int(acceptor)]

    if (donor_atom in ligand_atoms and acceptor_atom in protein_atoms) or (acceptor_atom in ligand_atoms and donor_atom in protein_atoms):
      if debug_info:
        # Additional information, print of the arrays obtained from MDAnalysis.
        print(
          f'Frame {int(frame)} | '
          f'{donor_atom.resname} -- {donor_atom.name} | '
          f'{acceptor_atom.resname} -- {acceptor_atom.name} | '
          f'd={dist:.2f} Å angle={angle:.1f}'
        )
      host_guest_hbonds.append([frame, donor, hydrogen, acceptor, dist, angle])

  print(f'\n - Keeping {len(host_guest_hbonds)}/{len(hbonds)} H-bonds (total amount, not unique H-bonds).')
  return host_guest_hbonds, protein_atoms

def hbond_population_counter(universe_mda,host_guest_hbonds):
  hbond_population = {}

  # Initial state of dictionary #
  frame, donor_idx, hydrogen_idx, acceptor_idx, dist, angle = host_guest_hbonds[0]

  donor_atom = universe_mda.atoms[int(donor_idx)]
  acceptor_atom = universe_mda.atoms[int(acceptor_idx)]
  hbond_population[f'hb{0}'] = {
        'donor' : f'{donor_atom.resname}{donor_atom.resid}-{donor_atom.name}',
        'donor_idx' : int(donor_idx),
        'acceptor' : f'{acceptor_atom.resname}{acceptor_atom.resid}-{acceptor_atom.name}',
        'acceptor_idx' : int(acceptor_idx),
        'frames' : 0
      }

  # Generating dictionary and counting the number of frames each H-bond is present #

  hi = 1
  for hb in host_guest_hbonds:
    frame, donor_idx, hydrogen_idx, acceptor_idx, dist, angle = hb

    donor_atom = universe_mda.atoms[int(donor_idx)]
    name_donor = f'{donor_atom.resname}{donor_atom.resid}-{donor_atom.name}'
    acceptor_atom = universe_mda.atoms[int(acceptor_idx)]
    name_acceptor = f'{acceptor_atom.resname}{acceptor_atom.resid}-{acceptor_atom.name}'

    found = False
    for bond in hbond_population:
      if name_donor == hbond_population[bond]['donor'] and name_acceptor == hbond_population[bond]['acceptor']:
        hbond_population[bond]['frames'] += 1
        found = True
        break
    if not found:
      hbond_population[f'hb{hi}'] = {
        'donor' : f'{donor_atom.resname}{donor_atom.resid}-{donor_atom.name}',
        'donor_idx' : int(donor_idx),
        'acceptor' : f'{acceptor_atom.resname}{acceptor_atom.resid}-{acceptor_atom.name}',
        'acceptor_idx' : int(acceptor_idx),
        'frames' : 1
      }
      hi += 1

  return hbond_population

def hbond_filter2(n_frames,step_hbond,universe_mda,population_hbond,host_guest_hbonds,debug_info):
  #* H-bond Filter 2: Population of H-bonds *#

  print(f'\nFilter 2: The H-bond must be in {population_hbond*100} % of the simulation.')

  hbond_population = hbond_population_counter(universe_mda,host_guest_hbonds)
  
  if debug_info:
    # Print to see the dictionary of the full unique hbonds before filtering by population.
    print('\nAll unique H-bonds:')
    console.print(Pretty(hbond_population))

  # Filtering H-bonds without enough population #

  keys_to_rm = [hbond for hbond in hbond_population if hbond_population[hbond]['frames'] <= (n_frames/step_hbond)*population_hbond]

  print(f' - Removing {len(keys_to_rm)}/{len(hbond_population)} H-bonds, they don\'t reach {population_hbond*100} % of the simulation..')

  for key in keys_to_rm:
    hbond_population.pop(key, 'None')

  if debug_info:
    # Print to see the dictionary of the full unique hbonds before filtering by population.
    print('\nH-bonds with more than half the simulation population:')
    console.print(Pretty(hbond_population))

  return hbond_population

def hbond_filter3(universe_mda,hbond_population,debug_info):
  #* H-bond Filter 3:  Filtering H-bonds of C-H*#

  print(f'\nFilter 3: Removing H-bonds where the acceptor or donor atom is a C.')
  keys_to_rm = []

  for hb in hbond_population:
    print(f'\n - {hb}')
    donor_idx = hbond_population[hb]['donor_idx']
    acceptor_idx = hbond_population[hb]['acceptor_idx']

    donor = universe_mda.select_atoms(f'index {donor_idx}')[0]
    print('   - Donor is:',donor.element)
    acceptor = universe_mda.select_atoms(f'index {acceptor_idx}')[0]
    print('   - Acceptor is:',acceptor.element)

    removed = False # to break the loop for the H-bond if it is already removed for any reason.
    if donor.element == 'H':
      neighbours = universe_mda.select_atoms(f"(around 2.0 index {donor.index} )and (resname {donor.resname} and resid {donor.resid})")# and protein")

      for n in neighbours:
        if 'C' in n.name:
          print(f'   - Removing {hb} because of C-H bond.')
          keys_to_rm.append(hb)
          removed = True
          break
        else:
          hbond_population[hb]['donor_idx'] = int(n.index)
          hbond_population[hb]['donor'] = f'{n.resname}{n.resid}-{n.name}'
          break

    elif donor.element == 'C' and not removed:
      print(f'   - Removing {hb} because donor is a C.')
      keys_to_rm.append(hb)
      removed = True

    if acceptor.element == 'H' and not removed:
      neighbours = universe_mda.select_atoms(f"(around 2.0 index {acceptor.index} )and (resname {acceptor.resname} and resid {acceptor.resid})")# and protein")

      for n in neighbours:
        if 'C' in n.name:
          print(f'   - Removing {hb} because of C-H bond.')
          keys_to_rm.append(hb)
          removed = True
          break
        else:
          hbond_population[hb]['acceptor_idx'] = int(n.index)
          hbond_population[hb]['acceptor'] = f'{n.resname}{n.resid}-{n.name}'
          break

    elif acceptor.element == 'C' and not removed:
      print(f'   - Removing {hb} because acceptor is a C.')
      keys_to_rm.append(hb)

  print(f'\n - Removing {len(keys_to_rm)}/{len(hbond_population)} H-bonds, C is the acceptor or donor atom.')

  for key in keys_to_rm:
    hbond_population.pop(key, 'None')

  if debug_info:
    # Print to see the dictionary of the full unique hbonds before filtering by population.
    print('\nH-bonds without a C as an acceptor or donor atom:')
    console.print(Pretty(hbond_population))

  return hbond_population

def hbond_search(universe_mda, guest_resname, ligand_atoms, step_hbond, population_hbond, d_DH_cutoff, d_AH_cutoff, debug_info):
  #! Step 2: H-bond !#

  #* Searching for H-bonds in the simulation *#

  n_frames = len(universe_mda.trajectory)

  print(f'\n 2. Searching for H-bonds in the simulation in {int(n_frames/step_hbond)} / {n_frames} frames.\n')

  hbond_analysis = HydrogenBondAnalysis(
    universe      = universe_mda,
    donors_sel    = f'protein or resname {guest_resname}',
    acceptors_sel = f'protein or resname {guest_resname}',
    hydrogens_sel = 'name H*',
    d_h_cutoff    = d_DH_cutoff,
    d_a_cutoff    = d_AH_cutoff
  )

  hbond_analysis.run(verbose=True, step=step_hbond)
  hbonds = hbond_analysis.results.hbonds

  host_guest_hbonds, protein_atoms = hbond_filter1(universe_mda, ligand_atoms, hbonds, debug_info)

  hbond_population = hbond_filter2(n_frames,step_hbond,universe_mda,population_hbond,host_guest_hbonds,debug_info)

  hbond_population = hbond_filter3(universe_mda,hbond_population,debug_info)

  return hbond_population, protein_atoms

def anchor_finder(mda_guest_candidates_idx,universe_mda,g0,hx):
  print('\n   - Selecting Triads of anchor candidates for Guest.')
  
  bonded_g0 = []
  bondeds_g1 = []
  triads_guest_atoms = []
  triads_host_atoms = [] 

  if g0.index+2 in mda_guest_candidates_idx:
    bonds_acceptor = g0.bonds
    for bond in bonds_acceptor:
      bonded_atom = bond.partner(g0)
      if bonded_atom.index+2 in mda_guest_candidates_idx:
        bonded_g0.append(bonded_atom)
        
    print(f'\n     - Bonded atoms to G0 ({g0.resname}{g0.resid}-{g0.name}) : {[atom.name for atom in bonded_g0]}')
    for g1 in bonded_g0:
      bonds_g1 = g1.bonds
      g1_temp = []
      for bond in bonds_g1:
        bonded_atom = bond.partner(g1)
        if bonded_atom != g0 and bonded_atom.index+2 in mda_guest_candidates_idx:
          g1_temp.append(bonded_atom)
      bondeds_g1.append(g1_temp)
      print(f'       - Bonded atoms to G1 ({g1.name}) : {[atom.name for atom in g1_temp]}')

    for i in range(len(bonded_g0)):
      g1 = bonded_g0[i]
      bonded_to_g1 = bondeds_g1[i]
      print(f'\n     - Found Triads of guest atoms for G0 ({g0.name}) and G1 ({g1.name})')
      for j in range(len(bonded_to_g1)):
        triad = [int(g0.index), int(g1.index), int(bonded_to_g1[j].index)]
        triads_guest_atoms.append(triad)
        print(f'       - Triad {j} : {[idx for idx in triad]}')
        print(f'                   {[universe_mda.atoms[idx].name for idx in triad]}')

    print('\n   - Selecting Triads of anchor candidates for Host.')
    print(f'     - H-bond in residue {hx.resname}{hx.resid}, selecting their CA, C and N.')

    alpha_c = universe_mda.select_atoms(f'name CA and resname {hx.resname} and resid {hx.resid}')[0]
    n = universe_mda.select_atoms(f'name N and resname {hx.resname} and resid {hx.resid}')[0]
    c = universe_mda.select_atoms(f'name C and resname {hx.resname} and resid {hx.resid}')[0]

    triads_host_atoms = [[int(n.index), int(alpha_c.index), int(c.index)],[int(c.index), int(alpha_c.index), int(n.index)]]
    print(f'     - Triads for the Host ({hx.resname}{hx.resid}) :')
    for i in range(len(triads_host_atoms)):
      print(f'       - Triad {i} : {[idx for idx in triads_host_atoms[i]]}')
      print(f'                   {[universe_mda.atoms[idx].name for idx in triads_host_atoms[i]]}')

  print(triads_guest_atoms)
  print(triads_host_atoms)
  triads_guest_atoms = [t for t in triads_guest_atoms if len(t) == 3]
  triads_host_atoms = [t for t in triads_host_atoms if len(t) == 3]

  return triads_guest_atoms, triads_host_atoms


def find_triads(universe_mda,hbond_population,mda_guest_candidates_idx,protein_atoms,debug_info):
  #! Step 3: Finding triads of atoms for guest and host !#

  print('\n 3. Searching anchor host-guest cominations. \n ')

  #* Creating list of related atoms*#

  print('\nStarting to look for restraints\'s anchor candidates.')

  candidates_restrains = []
  for hb in hbond_population:

    print(f'\n - Looking anchors using H-bond "{hb}".')

    donor = universe_mda.select_atoms(f'index {hbond_population[hb]['donor_idx']}')[0]
    acceptor = universe_mda.select_atoms(f'index {hbond_population[hb]['acceptor_idx']}')[0]

    if donor in protein_atoms:
      g0 = acceptor
      if g0.index+2 not in mda_guest_candidates_idx:
        g0 = g0.bonds[0].partner(g0)
      hx = donor
      
    else:
      g0 = donor
      if g0.index+2 not in mda_guest_candidates_idx:
        g0 = g0.bonds[0].partner(g0)
      hx = acceptor

    triads_guest_atoms, triads_host_atoms = anchor_finder(mda_guest_candidates_idx,universe_mda,g0,hx)

    print(f'\n   - Combining {len(triads_host_atoms)} triads of possible host anchors with {len(triads_guest_atoms)} triads of possible guest anchors. ')

    for i in range(len(triads_host_atoms)):
      for j in range(len(triads_guest_atoms)):
        candidates_restrains.append(triads_host_atoms[i]+triads_guest_atoms[j])

  print('\nCombining all sets of possible anchors.')

  if debug_info:
    print('All possible candidate sets:')
    # Additional information: We print all the candidates for restrains, there can be repeated ones.   
    console.print(Pretty(candidates_restrains))

  print('\nCleaning repeated sets of anchors.')

  unique_candidates_restraints = []
  for i in range(len(candidates_restrains)):
    if candidates_restrains[i] not in unique_candidates_restraints:
         unique_candidates_restraints.append(candidates_restrains[i])

  print('\nUnique sets of candidates for restrain anchors:')
  console.print(Pretty(unique_candidates_restraints))

  return unique_candidates_restraints

def compute_angles_restr(universe_mda,unique_candidates_restraints,set_idx):

  h0 = universe_mda.select_atoms(f' index {unique_candidates_restraints[set_idx][0]}')
  h1 = universe_mda.select_atoms(f' index {unique_candidates_restraints[set_idx][1]}')
  g0 = universe_mda.select_atoms(f' index {unique_candidates_restraints[set_idx][3]}')
  g1 = universe_mda.select_atoms(f' index {unique_candidates_restraints[set_idx][4]}')

  pos_h0 = h0.positions[0]
  pos_h1 = h1.positions[0]
  pos_g0 = g0.positions[0]
  pos_g1 = g1.positions[0]

  theta_A = float(np.degrees(calc_angles(pos_h1,pos_h0,pos_g0)))
  theta_B = float(np.degrees(calc_angles(pos_h0,pos_g0,pos_g1)))

  return theta_A, theta_B

def check_angles(universe_mda,min_angle,max_angle,unique_candidates_restraints,debug_info):
  #! Step 4: Checking angles !#

  print('\n 4. Filtering possible anchors by their angle. \n ')
  print('\n - We will compute theta_A (H1-H0-G0) and theta_B (H0-G0-G1) and check if they are')
  print(f'   in the range [{min_angle}º-{max_angle}º]\n')

  idx_to_remove = []
  for i in range(len(unique_candidates_restraints)):
    print(f' - Working on set {i}  {unique_candidates_restraints[i]}')
    print(f'                     {[universe_mda.atoms[idx].name for idx in unique_candidates_restraints[i]]}')

    theta_A, theta_B = compute_angles_restr(universe_mda,unique_candidates_restraints,i)

    removed = False

    print('\n   - Checking angles of the first frame.')
    if not (45 <= theta_A <= 135) or not (45 <= theta_B <= 135):
      print('     - Not keeping this set of candidates. First frame\'s theta_A or theta_B is out of [45º,135º].')
      print(f'       theta_A = {theta_A} º // theta_B = {theta_B} º\n')

      idx_to_remove.append(i)
      removed = True

    else:
      print('   - Checking angles for the rest of the simulation.')

      if debug_info:
        print('frame     theta_A     theta_B')

      for ts in universe_mda.trajectory:
        if removed:
          break

        theta_A, theta_B = compute_angles_restr(universe_mda,unique_candidates_restraints,i)

        if debug_info:
          print(ts.frame, theta_A, theta_B)

        if not (45 <= theta_A <= 135) or not (45 <= theta_B <= 135):
          print('     - Not keeping this set of candidates. theta_A or theta_B is out of [45º,135º] during he simulation.')
          print(f'       theta_A = {theta_A} º // theta_B = {theta_B} º \n')
          idx_to_remove.append(i)
          removed = True
    if not removed:
      print('   - Keeping this set of candidates\n')

  print(f' - Removing {len(idx_to_remove)}/{len(unique_candidates_restraints)} sets of candidates, angles out of range.')

  for idx in sorted(idx_to_remove, reverse=True):
    unique_candidates_restraints.pop(idx)

  if debug_info:
    print('\nSets of restraints after filtering angles out of range:')
    console.print(Pretty(unique_candidates_restraints))
  
  return unique_candidates_restraints

def check_distance_guest_COM(universe_mda,unique_candidates_restraints,ligand_atoms,debug_info):
  #! Step 5:  Selecting the ones with the closest G0 to the COM of the ligand !#

  print('\n 5. Keeping sets with G0 closest to Center of Mass of the Guest. \n')

  avg_dists = []
  print(' - Computing distance from G0 to guest\'s COM for every frame of simulation.')
  for i in range(len(unique_candidates_restraints)):
    g0 = universe_mda.select_atoms(f' index {unique_candidates_restraints[i][3]}')
    pos_g0 = g0.positions[0]

    dists = []
    for ts in universe_mda.trajectory:
      pos_g0 = g0.positions[0]
      com = ligand_atoms.center_of_mass()
      dist = np.linalg.norm(pos_g0 - com)
      dists.append(dist)
    if not dists:
        print(f"Skipping candidate {i}: no trajectory data")
        continue

    avg_dist = np.average(dists)

    print(f'   - G0 ({g0[0].resname}{g0[0].resid}-{g0[0].name}) distance to COM = {avg_dist.round(5)} Å')
    avg_dists.append(float(avg_dist))

  final_candidates = []
  minim = np.min(avg_dists)
  print(f'\nThe G0 with smallest distance to COM has = {minim.round(3)} Å')

  for i in range(len(avg_dists)):
    if round(avg_dists[i], 3) == round(np.min(avg_dists), 3):
      final_candidates.append(unique_candidates_restraints[i])

  print(f'\nWe are keeping {len(final_candidates)}/{len(unique_candidates_restraints)} sets of candidates.')

  if debug_info:
    # Additional Info: final sets of candidates.
    print(f'\nFinal sets of candidates:')
    console.print(Pretty(final_candidates))
  return final_candidates

def circ_std(angles):
  '''
  We use scipy.stats circvar to compute the circular standard deviation from the inputed angles in radians.
  '''
  V = circvar(angles, high=np.pi, low=-np.pi)
  R = 1 - V
  R = max(R, 1e-8) # if R goes to 0, then std would go to inf...
  return np.sqrt(-2 * np.log(R))

def compute_restraints_statistics(universe_mda,final_candidates,set_idx):
  print(f'   - Working on set {set_idx}  {final_candidates[set_idx]}')
  print(f'                       {[universe_mda.atoms[idx].name for idx in final_candidates[set_idx]]}')
  h0 = universe_mda.select_atoms(f' index {final_candidates[set_idx][0]}')
  h1 = universe_mda.select_atoms(f' index {final_candidates[set_idx][1]}')
  h2 = universe_mda.select_atoms(f' index {final_candidates[set_idx][2]}')
  g0 = universe_mda.select_atoms(f' index {final_candidates[set_idx][3]}')
  g1 = universe_mda.select_atoms(f' index {final_candidates[set_idx][4]}')
  g2 = universe_mda.select_atoms(f' index {final_candidates[set_idx][5]}')

  r_aA_time = []
  theta_A_time = []
  theta_B_time = []
  phi_A_time = []
  phi_B_time = []
  phi_C_time = []

  for ts in universe_mda.trajectory:
    pos_h0 = h0.positions[0]
    pos_h1 = h1.positions[0]
    pos_h2 = h2.positions[0]
    pos_g0 = g0.positions[0]
    pos_g1 = g1.positions[0]
    pos_g2 = g2.positions[0]

    r_aA = calc_bonds(pos_h0,pos_g0)
    r_aA_time.append(r_aA)

    theta_A = calc_angles(pos_h1,pos_h0,pos_g0)
    theta_A_time.append(theta_A)

    theta_B = calc_angles(pos_h0,pos_g0,pos_g1)
    theta_B_time.append(theta_B)

    phi_A = calc_dihedrals(pos_h2,pos_h1,pos_h0,pos_g0)
    phi_A_time.append(phi_A)

    phi_B = calc_dihedrals(pos_h1,pos_h0,pos_g0,pos_g1)
    phi_B_time.append(phi_B)

    phi_C = calc_dihedrals(pos_h0,pos_g0,pos_g1,pos_g2)
    phi_C_time.append(phi_C)

  avg_r_aA = float(np.average(r_aA_time))
  std_r_aA = float(np.std(r_aA_time))

  avg_theta_A = float(np.average(theta_A_time))
  std_theta_A = float(np.std(theta_A_time))

  avg_theta_B = float(np.average(theta_B_time))
  std_theta_B = float(np.std(theta_B_time))

  avg_phi_A = float(np.average(phi_A_time))
  std_phi_A = float(circ_std(phi_A_time))

  avg_phi_B = float(np.average(phi_B_time))
  std_phi_B = float(circ_std(phi_B_time))

  avg_phi_C = float(np.average(phi_C_time))
  std_phi_C = float(circ_std(phi_C_time))

  last_frame_vars = [float(r_aA),float(theta_A)*180/np.pi,float(theta_B)*180/np.pi,float(phi_A)*180/np.pi, float(phi_B)*180/np.pi,float(phi_C)*180/np.pi]

  return avg_r_aA, std_r_aA, avg_theta_A, std_theta_A, avg_theta_B, std_theta_B, avg_phi_A, std_phi_A, avg_phi_B, std_phi_B, avg_phi_C, std_phi_C,last_frame_vars

def scoring_candidates(universe_mda,final_candidates):

  #! Step 6 : Scoring the remaining candidates !#

  print('\n 6. Scoring the final sets of candidates, and choosing the anchors as the set with lowest score. \n')
  
  print(f' - Computing the parameters of the six Boresch restraints over the simulation.')
  scores = [] 
  last_frames_vars = []
  for i in range(len(final_candidates)):
    
    avg_r_aA, std_r_aA, avg_theta_A, std_theta_A, avg_theta_B, std_theta_B, avg_phi_A, std_phi_A, avg_phi_B, std_phi_B, avg_phi_C, std_phi_C, last_frame_vars = compute_restraints_statistics(universe_mda,final_candidates,i)
    last_frames_vars.append(last_frame_vars)

    score = (std_r_aA*std_theta_A*std_theta_B*std_phi_A*std_phi_B*std_phi_C*(avg_r_aA**2))/(np.sin(avg_theta_A)*np.sin(avg_theta_B))
    scores.append(float(score))

    print('\n                       Average          Standard deviation / circular st. dev.')
    print(f'     r_aA    :          {round(avg_r_aA,5)}          {round(std_r_aA,5)}')
    print(f'     theta_A :          {round(avg_theta_A,5)}          {round(std_theta_A,5)}')
    print(f'     theta_B :          {round(avg_theta_B,5)}          {round(std_theta_B,5)}')
    print(f'     phi_A   :          {round(avg_phi_A,5)}          {round(std_phi_A,5)}')
    print(f'     phi_B   :          {round(avg_phi_B,5)}          {round(std_phi_B,5)}')
    print(f'     phi_C   :          {round(avg_phi_C,5)}          {round(std_phi_C,5)}')

    print(f'\n   - Score of set {i}  {score}\n')

  lowest_score_idx = scores.index(np.min(scores))

  print(f'\n\n - We are choosing as restraints set {lowest_score_idx}, score = {round(scores[lowest_score_idx],6)}')
  print(' - The anchors are :\n')
  atom_roles = ['H0','H1','H2','G0','G1','G2']
  for i in range(len(final_candidates[lowest_score_idx])):
    atom = universe_mda.atoms[final_candidates[lowest_score_idx][i]]
    print(f'    {atom_roles[i]} ---> {atom.resname}{atom.resid}-{atom.name} (index {atom.index+1})')

  return final_candidates[lowest_score_idx], last_frames_vars[lowest_score_idx]

def restraint_search_william(guest_sdf_name:str,
                             debug_info:bool,
                             pdb_name_search:str,
                             traj_name:str,
                             guest_resname:str,
                             step_hbond:int = 1, 
                             population_hbond:float = 0.5,
                             d_DH_cutoff:float = 1.35, 
                             d_AH_cutoff:float = 3.3,
                             min_angle:float = 45,
                             max_angle:float = 135
                             ):
  """Full procedure
  """
  guest_candidates_ids = find_guest_candidates(guest_sdf_name,debug_info)

  mda_guest_candidates_idx, universe_mda, ligand_atoms = setting_up_mda(pdb_name_search,traj_name,guest_resname,guest_candidates_ids)

  hbond_population, protein_atoms = hbond_search(universe_mda, guest_resname, ligand_atoms, step_hbond, population_hbond, d_DH_cutoff, d_AH_cutoff, debug_info)

  unique_candidates_restraints = find_triads(universe_mda,hbond_population,mda_guest_candidates_idx,protein_atoms,debug_info)

  unique_candidates_restraints = check_angles(universe_mda,min_angle,max_angle,unique_candidates_restraints,debug_info)

  final_candidates = check_distance_guest_COM(universe_mda,unique_candidates_restraints,ligand_atoms,debug_info)

  anchors,last_frame_vars = scoring_candidates(universe_mda,final_candidates)
  return anchors,universe_mda,last_frame_vars

#? main runner

def main():
  #* System and file names *#
  guest_sdf_name = '../test_with_OMM/2u.sdf'
  guest_resname  = 'UNK'
  pdb_name       = '../test_with_OMM/p38_2u_SOLV.pdb'
  traj_name      = '../test_with_OMM/omm_simulations/rep_0/1ns_traj.dcd'

  #* Extra information needed? *#
  debug_info = False

  #* Pretty print *#
  console = Console(width=200)  # or 240, 300, etc.

  #* H-bond search parameters *#
  step_hbond       = 10
  d_DH_cutoff      = 1.35 # Maximum distance (angstrom) from the Donor to the H. D-H --- A
  d_AH_cutoff      = 3.3  # Maximum distance (angstrom) from the Acceptor to the H. D-H --- A
  population_hbond = 0.5  # Percentage of the simulation where the hbond has to be present.


  #* Angle filter *#
  min_angle = 45  # in degrees
  max_angle = 135 # in degrees

  guest_candidates_ids = find_guest_candidates(guest_sdf_name,debug_info)

  mda_guest_candidates_idx, universe_mda, ligand_atoms = setting_up_mda(pdb_name,traj_name,guest_resname,guest_candidates_ids)

  hbond_population, protein_atoms = hbond_search(universe_mda, guest_resname, ligand_atoms, step_hbond, population_hbond, d_DH_cutoff, d_AH_cutoff, debug_info)

  unique_candidates_restraints = find_triads(universe_mda,hbond_population,mda_guest_candidates_idx,protein_atoms,debug_info)

  unique_candidates_restraints = check_angles(universe_mda,min_angle,max_angle,unique_candidates_restraints,debug_info)

  final_candidates = check_distance_guest_COM(universe_mda,unique_candidates_restraints,ligand_atoms,debug_info)

  anchors = scoring_candidates(universe_mda,final_candidates)

#? execution

if __name__ == '__main__':
  main()