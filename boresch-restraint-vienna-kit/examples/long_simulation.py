from boresch_restraint_vienna_kit.simulation_julian import *
from boresch_restraint_vienna_kit.restraints_william import *
from boresch_restraint_vienna_kit.restraints_openfe import *
from boresch_restraint_vienna_kit.drawing_boresch_restraints import *
import MDAnalysis as mda

import openmm.unit as unit

pdb_name='2QBS'
sdf_orig_name='2qbs_B_024'
lig_resname='024'

preparing_from_PDBank(pdb_name,
                      lig_resname=lig_resname,
                      sdf_orig_name = sdf_orig_name,
                      residues_to_remove = {'CL', 'NA', lig_resname},
                      protonate=True)

supp = Chem.SDMolSupplier(f'lig_{lig_resname}.sdf', removeHs=False)
mol = supp[0]

def deprotonate_carboxylic_acids(mol):
    rw = Chem.RWMol(mol)
    
    # SMARTS: carboxylic acid
    pattern = Chem.MolFromSmarts("C(=O)[OH]")
    matches = mol.GetSubstructMatches(pattern)

    for match in matches:
        carbonyl_c = match[0]
        double_bond_o = match[1]
        hydroxyl_o = match[2]

        o_atom = rw.GetAtomWithIdx(hydroxyl_o)

        # remove hydrogen(s) on OH oxygen (if explicit)
        for nbr in list(o_atom.GetNeighbors()):
            if nbr.GetAtomicNum() == 1:  # hydrogen
                rw.RemoveAtom(nbr.GetIdx())

        # set formal charge on oxygen
        o_atom.SetFormalCharge(-1)

    return rw.GetMol()

deprot_mol = deprotonate_carboxylic_acids(mol)
Chem.SanitizeMol(deprot_mol)


w = Chem.SDWriter(f'lig_{lig_resname}_deprot.sdf')
w.write(deprot_mol)
w.close()


preparation_simulations(molecule_name            = f'lig_{lig_resname}_deprot' ,
                        pdb_name                 = f'{pdb_name}_CLEAN.pdb',
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
                        traj_numb         = 3,
                        md_prod_total     = 150.0*unit.nanosecond,
                        ini               = 1,
                        fin               = 3
                        )