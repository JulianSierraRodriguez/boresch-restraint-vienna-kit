from boresch_restraint_vienna_kit.simulation_julian import *
from boresch_restraint_vienna_kit.restraints_william import *
from boresch_restraint_vienna_kit.restraints_openfe import *
from boresch_restraint_vienna_kit.drawing_boresch_restraints import *
import MDAnalysis as mda
import openmm.unit as unit
import pathlib, os
import pandas as pd

number_repeats = 10


host_anchor_list_W = []
guest_anchor_list_W = []
host_names_list_W = []
guest_names_list_W = []
protein_res_W = []
r_aA_W = []
theta_A_W = []
theta_B_W = []
phi_A_W = []
phi_B_W = []
phi_C_W = []

host_anchor_list_openfe = []
guest_anchor_list_openfe = []
host_names_list_openfe = []
guest_names_list_openfe = []
protein_res_openfe = []
r_aA_openfe = []
theta_A_openfe = []
theta_B_openfe = []
phi_A_openfe = []
phi_B_openfe = []
phi_C_openfe = []

for i in range(number_repeats):
  folder = f'rep_{i}'
  path = pathlib.Path(folder)
  path.mkdir(exist_ok=True)
  os.system(f'cp *.pdb {folder}')
  os.system(f'cp *.sdf {folder}')
  os.chdir(folder)
  


  preparation_simulations(molecule_name            = '2f' ,
                          pdb_name                 = 'p38_protein_CLEAN.pdb',
                          solvated_PDB             = False,
                          folder_prep              = 'system_prep',
                          timestep                 = 2.0*unit.femtosecond,
                          friction_coeff           = 1.0*1/unit.picosecond,
                          md_temperature           = 298.15*unit.kelvin,
                          minimization_max_steps   = 5000,
                          min_tolerance            = 4 * unit.kilojoule_per_mole / unit.nanometer,
                          pdb_minimization_final   = 'minimized_system.pdb',
                          reporter_interval        = 5_000,
                          temp_increment           = 100, # mine 5 or openfe (no ramp) 298.15  
                          md_equil_nvt_length      = 0.25*unit.nanosecond,
                          md_steps_after_ramp      = 0,      #mine 5000 or openfe 0
                          pdb_NVT_final            = 'NVT_system.pdb',
                          md_equil_npt_length      = 0.5*unit.nanosecond,
                          md_pressure              = 1*unit.bar,
                          pdb_NPT_final            = 'NPT_system.pdb'
                          )

  production_simulations( folder_prep       = 'system_prep',
                          timestep          = 2.0*unit.femtosecond,
                          friction_coeff    = 1.0*1/unit.picosecond,
                          md_temperature    = 298.15*unit.kelvin,
                          reporter_interval = 5_000,
                          traj_interval     = 1_000,
                          traj_numb         = 1,
                          md_prod_total     = 1.0*unit.nanosecond,
                          ini               = 1,
                          fin               = 1
                          )

  anchors, universe_mda, last_frame_vars = restraint_search_william(guest_sdf_name   = '2f.sdf',
                                                   debug_info       = False,
                                                   pdb_name_search  = 'last_frame_1.pdb',
                                                   traj_name        = 'traj_1ns_1.dcd',
                                                   guest_resname    = 'UNK',
                                                   step_hbond       = 1, 
                                                   population_hbond = 0.5,
                                                   d_DH_cutoff      = 1.35, 
                                                   d_AH_cutoff      = 3.3,
                                                   min_angle        = 45,
                                                   max_angle        = 135
                                                   )

  resname_L, resname_P, resid_L, resid_P, atom_names = williams_anchors_to_names(universe_mda, anchors)

  host_anchor_list_W.append(anchors[:3])
  guest_anchor_list_W.append(anchors[3:])
  host_names_list_W.append(atom_names[:3])
  guest_names_list_W.append(atom_names[3:])
  protein_res_W.append(f'{resname_P}{resid_P}')
  r_aA_W.append(last_frame_vars[0])
  theta_A_W.append(last_frame_vars[1])
  theta_B_W.append(last_frame_vars[2])
  phi_A_W.append(last_frame_vars[3])
  phi_B_W.append(last_frame_vars[4])
  phi_C_W.append(last_frame_vars[5])

  

  drawing(path_pdb      = 'last_frame_1.pdb', 
          path_lig_sdf  = '2f.sdf', 
          resname_lig   = resname_L,
          resid_lig     = resid_L, 
          resname_prot  = resname_P,
          resid_prot    = resid_P,
          anchors_names = atom_names,
          figure_name   = f'complex_william_{folder}'
              )

  #we dont want to use this universe because the ligand is no longer 1, openfe makes it 351 (protein has 350)
  u, host_atoms, guest_atoms, last_frame_vars = restraint_search_openfe(pdb_name_search = 'last_frame_1.pdb',
                                                           guest_sdf_name  = '2f.sdf',
                                                           guest_resname   = 'UNK',
                                                           trajectory_name = 'traj_1ns_1.dcd',
                                                           path            = '.',
                                                           temperature     = 298.15 
                                                           ) 
  

  resname_L, resname_P, resid_L, resid_P, atom_names = openfe_anchors_to_names(universe_mda, host_atoms,guest_atoms)

  host_anchor_list_openfe.append(host_atoms)
  guest_anchor_list_openfe.append(guest_atoms)
  host_names_list_openfe.append(atom_names[:3])
  guest_names_list_openfe.append(atom_names[3:])
  protein_res_openfe.append(f'{resname_P}{resid_P}')
  r_aA_openfe.append(last_frame_vars[0])
  theta_A_openfe.append(last_frame_vars[1])
  theta_B_openfe.append(last_frame_vars[2])
  phi_A_openfe.append(last_frame_vars[3])
  phi_B_openfe.append(last_frame_vars[4])
  phi_C_openfe.append(last_frame_vars[5])

  drawing(path_pdb      = 'last_frame_1.pdb', 
          path_lig_sdf  = '2f.sdf', 
          resname_lig   = resname_L,
          resid_lig     = resid_L, 
          resname_prot  = resname_P,
          resid_prot    = resid_P,
          anchors_names = atom_names,
          figure_name   = f'complex_openfe_{folder}'
          )
  
  os.chdir('..')

df_openfe = pd.DataFrame({
    "method": "openfe",
    "host_anchor": host_anchor_list_openfe,
    "guest_anchor": guest_anchor_list_openfe,
    "host_names": host_names_list_openfe,
    "guest_names": guest_names_list_openfe,
    "protein_res": protein_res_openfe,
    "r_aA": r_aA_openfe,
    "theta_A": theta_A_openfe,
    "theta_B": theta_B_openfe,
    "phi_A": phi_A_openfe,
    "phi_B": phi_B_openfe,
    "phi_C": phi_C_openfe
})

df_W = pd.DataFrame({
    "method": "W",
    "host_anchor": host_anchor_list_W,
    "guest_anchor": guest_anchor_list_W,
    "host_names": host_names_list_W,
    "guest_names": guest_names_list_W,
    "protein_res": protein_res_W,
    "r_aA": r_aA_W,
    "theta_A": theta_A_W,
    "theta_B": theta_B_W,
    "phi_A": phi_A_W,
    "phi_B": phi_B_W,
    "phi_C": phi_C_W
})

df = pd.concat([df_openfe, df_W], ignore_index=True)

df.to_csv("comparison_W_openfe_10x1ns.csv", index=False)
df.to_excel("comparison_W_openfe_10x1ns.xlsx", index=False)
df.to_pickle("comparison_W_openfe_10x1ns.pkl")
