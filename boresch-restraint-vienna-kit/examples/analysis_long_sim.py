from boresch_restraint_vienna_kit.plot_restraints_mda import *
from boresch_restraint_vienna_kit.restraints_william import *
from boresch_restraint_vienna_kit.drawing_boresch_restraints import *
from boresch_restraint_vienna_kit.restraints_openfe import *
import MDAnalysis as mda

ref_pdb = 'system_prep/last_frame_0.pdb'
ref_sdf = 'system_prep/lig_024_deprot.sdf'

u = mda.Universe(
    ref_pdb,
    'traj_50ns_1.dcd'
  )

print(f' We have a dt = {u.trajectory.dt} ps')
dt = u.trajectory.dt
n_frames_1ns = int(1000/dt)

print(f'We have {n_frames_1ns} frames in 1ns')

with mda.Writer("first_1ns.dcd", n_atoms=u.atoms.n_atoms) as W:
  for ts in u.trajectory[:n_frames_1ns]:
    W.write(u.atoms)


anchors, universe_mda, last_frame_vars = restraint_search_william(guest_sdf_name   = ref_sdf,
                                                                  debug_info       = False,
                                                                  pdb_name_search  = ref_pdb,
                                                                  traj_name        = f'first_1ns.dcd',
                                                                  guest_resname    = 'UNK',
                                                                  step_hbond       = 1, 
                                                                  population_hbond = 0.5,
                                                                  d_DH_cutoff      = 1.35, 
                                                                  d_AH_cutoff      = 3.3,
                                                                  min_angle        = 45,
                                                                  max_angle        = 135
                                                   )

resname_L, resname_P, resid_L, resid_P, atom_names = williams_anchors_to_names(universe_mda, anchors)  

drawing(path_pdb      = ref_pdb, 
        path_lig_sdf  = ref_sdf, 
        resname_lig   = resname_L,
        resid_lig     = resid_L, 
        resname_prot  = resname_P,
        resid_prot    = resid_P,
        anchors_names = atom_names,
        figure_name   = f'complex_william'
              )

plot_restraints(pdb_name = ref_pdb,
                traj_name = ['traj_50ns_1.dcd','traj_50ns_2.dcd','traj_50ns_3.dcd'],
                step = 10,
                anchors = anchors,
                several_simulations = False,
                total_simulations = 1,
                host_idx_corrector = 1,
                guest_idx_corrector = 2,
                references=False,
                figure_name = 'ptp1b_024_W')

u, host_atoms, guest_atoms, last_frame_vars = restraint_search_openfe(pdb_name_search = 'last_frame_1.pdb',
                                                           guest_sdf_name  = '2f.sdf',
                                                           guest_resname   = 'UNK',
                                                           trajectory_name = 'first_1ns.dcd',
                                                           path            = '.',
                                                           temperature     = 298.15 
                                                           ) 
  

resname_L, resname_P, resid_L, resid_P, atom_names = openfe_anchors_to_names(universe_mda, host_atoms,guest_atoms)

anchors = host_atoms + guest_atoms

plot_restraints(pdb_name = ref_pdb,
                traj_name = ['traj_50ns_1.dcd','traj_50ns_2.dcd','traj_50ns_3.dcd'],
                step = 10,
                anchors = anchors,
                several_simulations = False,
                total_simulations = 1,
                host_idx_corrector = 1,
                guest_idx_corrector = 2,
                references=False,
                figure_name = 'ptp1b_024_openfe')


print('\n END of SCRIPT \n')
