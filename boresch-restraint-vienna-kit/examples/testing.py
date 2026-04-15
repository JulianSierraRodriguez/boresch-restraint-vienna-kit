from boresch_restraint_vienna_kit.simulation_julian import *
from boresch_restraint_vienna_kit.restraints_william import *
from boresch_restraint_vienna_kit.restraints_openfe import *
from boresch_restraint_vienna_kit.drawing_boresch_restraints import *
import MDAnalysis as mda

import openmm.unit as unit

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

traj_name = production_simulations( folder_prep       = 'system_prep',
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

print(last_frame_vars)

resname_L, resname_P, resid_L, resid_P, atom_names = williams_anchors_to_names(universe_mda, anchors)

drawing(path_pdb      = 'last_frame_1.pdb', 
        path_lig_sdf  = '2f.sdf', 
        resname_lig   = resname_L,
        resid_lig     = resid_L, 
        resname_prot  = resname_P,
        resid_prot    = resid_P,
        anchors_names = atom_names,
        figure_name   = 'complex_william'
            )

#we dont want to use this universe because the ligand is no longer 1, openfe makes it 351 (protein has 350)
u, host_atoms, guest_atoms, last_frame_vars = restraint_search_openfe(pdb_name_search = 'last_frame_1.pdb',
                                                         guest_sdf_name  = '2f.sdf',
                                                         guest_resname   = 'UNK',
                                                         trajectory_name = 'traj_1ns_1.dcd',
                                                         path            = '.',
                                                         temperature     = 298.15 
) 
print(last_frame_vars)


resname_L, resname_P, resid_L, resid_P, atom_names = openfe_anchors_to_names(universe_mda, host_atoms,guest_atoms)


drawing(path_pdb      = 'last_frame_1.pdb', 
        path_lig_sdf  = '2f.sdf', 
        resname_lig   = resname_L,
        resid_lig     = resid_L, 
        resname_prot  = resname_P,
        resid_prot    = resid_P,
        anchors_names = atom_names,
        figure_name   = 'complex_openfe'
        )

number_repeats = 1
lengthmd = traj_name.split('_')
print(f"comparison_W_openfe_{number_repeats}x{lengthmd}.csv")