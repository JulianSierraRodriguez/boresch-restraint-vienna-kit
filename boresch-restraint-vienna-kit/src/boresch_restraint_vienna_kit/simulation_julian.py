import openfe
from openmm import *
import openmm.app as app
import openmm.unit as unit
from openmm import XmlSerializer
from rdkit import Chem
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
# from rich.pretty import Pretty
# from rich.console import Console
from sys import stdout
import numpy as np
import os, pathlib
# from pathlib import Path

# from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals

# from openfe.protocols.openmm_afe.equil_afe_settings import BoreschRestraintSettings
# from openfe.protocols.openmm_afe.abfe_units import ABFEComplexSetupUnit
# from openff.units import unit as off_unit
from openmm.app import PDBFile
# from pint import UnitRegistry



def preparing_ligand(molecule_name:str):
  #! Getting Ligand !#
  supp = Chem.SDMolSupplier(f"{molecule_name}.sdf", removeHs=False)
  orig_rdkit = supp[0]
  ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]

  #! Charging the Ligand !#
  # To have the ligand's FF parameters.

  charge_settings = OpenFFPartialChargeSettings(partial_charge_method="am1bcc", off_toolkit_backend="ambertools")

  ligands = bulk_assign_partial_charges(
    molecules=ligands,
    overwrite=False,  
    method=charge_settings.partial_charge_method,
    toolkit_backend=charge_settings.off_toolkit_backend,
    generate_n_conformers=charge_settings.number_of_conformers,
    nagl_model=charge_settings.nagl_model,
    processors=1
  )
  return orig_rdkit,ligands

def system_solvated_PDB(pdb_name:str,
                        orig_rdkit,
                        forcefield:list = ["amber/ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml", "amber/phosaa10.xml"],
                        forcefield_small_molecule:str = "openff-2.2.1" 
                        ):
  #! We will use the coordinates from the PDB !#

    pdb = app.PDBFile(pdb_name)
    positions = pdb.positions
    topology = pdb.topology

    off_mol = Molecule.from_rdkit(orig_rdkit, allow_undefined_stereo=True)

    modeller = app.Modeller(topology, positions)

    #? We are using a PDB that was generated with OpenFE, the minimized.pdb, then we need to generate the Periodic Box. 

    # We obtain the positions of all atoms.
    coords = np.array([[x.value_in_unit(unit.nanometer),
                        y.value_in_unit(unit.nanometer),
                        z.value_in_unit(unit.nanometer)] for x,y,z in positions])

    # From all positions we get the minimum/maximum coordinates between all the positions.
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)

    # It is good to add some padding, if it is too much during the equilibration of density (NPT) the box will decrease its volume.

    box_lengths = max_coords - min_coords + 3.0  # add 2 nm padding

    modeller.topology.setUnitCellDimensions(unit.Quantity(box_lengths, unit.nanometer))
    modeller.positions = positions

    system_generator = SystemGenerator(
      forcefields=forcefield,
      small_molecule_forcefield=forcefield_small_molecule,
      molecules=[off_mol],
      periodic_forcefield_kwargs={
        'nonbondedMethod': app.PME,
        'nonbondedCutoff': 0.9*unit.nanometer,
        'constraints': app.HBonds,
        'rigidWater': True,
        'hydrogen_mass' :3.0 * unit.amu
      },
    )
    #* Adding the Periodic Box Vectors obtained previously to have periodicity. *#

    system = system_generator.create_system(modeller.topology)
    system.setDefaultPeriodicBoxVectors(
      Vec3(box_lengths[0], 0, 0),
      Vec3(0, box_lengths[1], 0),
      Vec3(0, 0, box_lengths[2])
    )
    return system, topology, positions

def solvating_system(pdb_name,
                     ligands,
                     forcefield:list = ["amber/ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml", "amber/phosaa10.xml"],
                     forcefield_small_molecule:str = "openff-2.2.1",
                     protein_ligand_pdb = "protein_ligand_system.pdb",
                     solvated_pdb = "solvated_system.pdb",
):
  protein = app.PDBFile(pdb_name)

  off_mol = ligands[0].to_openff()
  ligand_top = off_mol.to_topology().to_openmm()
  ligand_pos = off_mol.conformers[0].to_openmm()

  modeller = app.Modeller(protein.topology, protein.positions)
  modeller.add(ligand_top,ligand_pos)

  with open(protein_ligand_pdb, "w") as f:
    app.PDBFile.writeFile(
      modeller.topology,
      modeller.positions,
      f,
      keepIds=True
    )

  system_generator = SystemGenerator(
    forcefields=forcefield,
    small_molecule_forcefield=forcefield_small_molecule,
    molecules=[off_mol],
    periodic_forcefield_kwargs={
      'nonbondedMethod': app.PME,
      'nonbondedCutoff': 0.9*unit.nanometer,
      'constraints': app.HBonds,
      'rigidWater': True
    },
  )

  ff = system_generator.forcefield

  modeller.addSolvent(
    ff,
    model='tip3p',
    padding=1.0*unit.nanometer,
    ionicStrength = 0.15 * unit.molar
  )

  with open(solvated_pdb, "w") as f:
    app.PDBFile.writeFile(
      modeller.topology,
      modeller.positions,
      f,
      keepIds=True
    )
  system = system_generator.create_system(modeller.topology)
  topology = modeller.topology
  positions = modeller.positions
  return system, topology, positions

def setting_up_simulation(md_temperature: unit.Quantity,
                          friction_coeff:unit.Quantity,
                          timestep:unit.Quantity,
                          topology,
                          system,
                          positions):
  integrator = LangevinMiddleIntegrator(
    md_temperature,
    friction_coeff,
    timestep
  )

  platform = Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'mixed'}
  simulation = app.Simulation(
    topology,
    system,
    integrator,
    platform,
    properties
  )
  simulation.context.setPositions(positions)
  return integrator, simulation

def minimization(simulation,
                 min_tolerance,
                 minimization_max_steps,
                 pdb_minimization_final
                 ):
  potential_energy_ini = simulation.context.getState(getEnergy=True).getPotentialEnergy()
  print(f'\nInitial potential energy: {potential_energy_ini}')

  simulation.minimizeEnergy(
    tolerance = min_tolerance,
    maxIterations= minimization_max_steps
  )

  potential_energy_fin = simulation.context.getState(getEnergy=True).getPotentialEnergy()
  print(f'\nFinal potential energy: {potential_energy_fin}')

  with open(pdb_minimization_final, 'w') as f:
    app.PDBFile.writeFile(
      simulation.topology,
      simulation.context.getState(getPositions=True).getPositions(),
      file=f
    )
  return simulation

def simulation_NVT(simulation,
                   integrator,
                   reporter_interval:int, 
                   md_equil_nvt_length:unit.Quantity,
                   timestep:unit.Quantity,
                   md_temperature:unit.Quantity,
                   temp_increment:int,
                   md_steps_after_ramp:int,
                   pdb_heating_final:str
                   ):
  simulation.reporters.append(
    app.StateDataReporter(
      stdout,
      reporter_interval,
      step=True, potentialEnergy=True, temperature=True, volume=True, density=True
    )
  )

  print('\nHeating the system - Running NVT equilibration')
  simulation.context.setVelocitiesToTemperature(5*unit.kelvin)

  mdsteps = round(md_equil_nvt_length/timestep)
  if float(int((md_temperature/temp_increment)/unit.kelvin)) == round(((md_temperature/temp_increment)/unit.kelvin),1):
    stages = int((md_temperature/temp_increment)/unit.kelvin)
  else :
    stages = int((md_temperature/temp_increment)/unit.kelvin)+1
  steps_per_stage = int(mdsteps//stages)

  for i in range(stages):
    temperature = (temp_increment+(i*temp_increment))*unit.kelvin
    if i%20 == 0:
      print(f'Heating to {temperature}')
    integrator.setTemperature(temperature)
    simulation.step(steps_per_stage)

  # Now we should be at 298.15K, but we will add some more steps to maintain T equilibrated
  print(f'The temperature now should be around {temperature}.')
  simulation.step(md_steps_after_ramp)

  with open(pdb_heating_final, 'w') as f:
    app.PDBFile.writeFile(
      simulation.topology,
      simulation.context.getState(getPositions=True).getPositions(),
      file=f
    )
  return simulation

def simulation_NPT(simulation,
                   system,
                   timestep:unit.Quantity       = 2.0*unit.femtosecond,
                   friction_coeff:unit.Quantity = 1.0*1/unit.picosecond,
                   md_temperature:unit.Quantity = 298.15*unit.kelvin,
                   md_equil_npt_length:unit.Quantity = 0.5*unit.nanosecond,
                   md_pressure:unit.Quantity = 1*unit.bar,
                   pdb_density_equil_final:str = 'NPT_system.pdb'
                   ):
  integrator = LangevinMiddleIntegrator(
      md_temperature,
      friction_coeff,
      timestep
    )

  mdsteps = round(md_equil_npt_length/timestep)
  barostat = system.addForce(MonteCarloBarostat(md_pressure, md_temperature))
  simulation.context.reinitialize(True)
  print('\nRunning NPT equilibration')
  simulation.step(mdsteps)

  npt_state = simulation.context.getState(
      getPositions=True,
      getVelocities=True,
      getEnergy=True,
      enforcePeriodicBox=True
  )

  with open(pdb_density_equil_final, 'w') as f:
    app.PDBFile.writeFile(
      simulation.topology,
      simulation.context.getState(getPositions=True).getPositions(),
      file=f
    )
  return simulation,system,npt_state

def save_checkpoint_for_restart(iter,system,simulation,integrator):
  with open(f"system_{iter}.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))
  with open(f"integrator_{iter}.xml", "w") as f:
    f.write(XmlSerializer.serialize(integrator))
  state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=False,
    enforcePeriodicBox=True
  )
  with open(f"state_{iter}.xml", "w") as f:
    f.write(XmlSerializer.serialize(state))

  simulation.saveCheckpoint(f'checkpoint_{iter}.chk')
  positions = state.getPositions()

  with open(f"last_frame_{iter}.pdb", "w") as f:
      PDBFile.writeFile(simulation.topology, positions, f)


def preparation_simulations(molecule_name:str,
                            pdb_name:str,
                            solvated_PDB:bool,
                            folder_prep:str              = 'system_prep',
                            timestep:unit.Quantity       = 2.0*unit.femtosecond,
                            friction_coeff:unit.Quantity = 1.0*1/unit.picosecond,
                            md_temperature:unit.Quantity = 298.15*unit.kelvin,
                            minimization_max_steps       = 5000,
                            min_tolerance:unit.Quantity  = 4 * unit.kilojoule_per_mole / unit.nanometer,
                            pdb_minimization_final       = 'minimized_system.pdb',
                            reporter_interval:int        = 5_000,
                            temp_increment:int           = 100, # mine 5 or openfe (no ramp) 298.15  
                            md_equil_nvt_length:unit.Quantity  = 0.25*unit.nanosecond,
                            md_steps_after_ramp:unit.Quantity  = 0,      #mine 5000 or openfe 0
                            pdb_NVT_final:str    = 'NVT_system.pdb',
                            md_equil_npt_length:unit.Quantity = 0.5*unit.nanosecond,
                            md_pressure:unit.Quantity = 1*unit.bar,
                            pdb_NPT_final:str = 'NPT_system.pdb'

                            ):

  path = pathlib.Path(f'./{folder_prep}')
  path.mkdir(exist_ok=True)
  os.chdir(folder_prep)

  os.system(f'cp ../{pdb_name} .')
  os.system(f'cp ../{molecule_name}.sdf .')

  orig_rdkit,ligands = preparing_ligand(molecule_name)
  
  if solvated_PDB:
    system, topology, positions = system_solvated_PDB(pdb_name,
                                                    orig_rdkit,
                                                    forcefield = ["amber/ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml", "amber/phosaa10.xml"],
                                                    forcefield_small_molecule = "openff-2.2.1" )

  else:
    system, topology, positions = solvating_system(pdb_name,
                                                 ligands,
                                                 forcefield = ["amber/ff14SB.xml", "amber/tip3p_standard.xml", "amber/tip3p_HFE_multivalent.xml", "amber/phosaa10.xml"],
                                                 forcefield_small_molecule = "openff-2.2.1",
                                                 protein_ligand_pdb = "protein_ligand_system.pdb",
                                                 solvated_pdb = "solvated_system.pdb",)

  integrator, simulation = setting_up_simulation(md_temperature,
                          friction_coeff,
                          timestep,
                          topology,
                          system,
                          positions)

  #* Minimizing system *#

  simulation = minimization(simulation,
                 min_tolerance,
                 minimization_max_steps,
                 pdb_minimization_final
                 )


  #* NVT equilibrium --- Heating *#
  # Changing the temperature slowly and having steps after the change of temperature to let the system slowly equilibrate as it gets heated.

  simulation = simulation_NVT(simulation,
                              integrator,
                              reporter_interval, 
                              md_equil_nvt_length,
                              timestep,
                              md_temperature,
                              temp_increment,
                              md_steps_after_ramp,
                              pdb_NVT_final)

  #* NPT equilibrium --- Density *#
  # We introduce a barostat and we allow the system to change its volume to equilibrate the density to around one.
  
  simulation,system,npt_state = simulation_NPT(simulation, 
                                               system,
                                               timestep,
                                               friction_coeff,
                                               md_temperature,
                                               md_equil_npt_length,
                                               md_pressure,
                                               pdb_NPT_final)
  
  save_checkpoint_for_restart(0,system,simulation,integrator)

  os.chdir('..')

def set_up_with_checkpoint(iter,
                           folder_prep,
                           timestep:unit.Quantity       = 2.0*unit.femtosecond,
                           friction_coeff:unit.Quantity = 1.0*1/unit.picosecond,
                           md_temperature:unit.Quantity = 298.15*unit.kelvin,
):
  
  if iter == 1:
    folder = folder_prep
  else:
    folder = '.'
  
  pdb = app.PDBFile(f'{folder}/last_frame_{iter-1}.pdb')

  with open(f"{folder}/system_{iter-1}.xml") as f:
      system = XmlSerializer.deserialize(f.read())

  with open(f"{folder}/integrator_{iter-1}.xml") as f:
      integrator = XmlSerializer.deserialize(f.read())

  platform = Platform.getPlatformByName("CUDA")
  properties = {"Precision": "mixed"}  # MUST match original run
  simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

  with open(f"{folder}/checkpoint_{iter-1}.chk", "rb") as f:
      simulation.context.loadCheckpoint(f.read())

  # Doing this seems to be stabilizing CUDA
  state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    enforcePeriodicBox=True
  )

  integrator = LangevinMiddleIntegrator(
    md_temperature,
    friction_coeff,
    timestep
  )
  platform = Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'mixed'}

  topology = pdb.topology
  simulation = app.Simulation(
    topology,
    system,
    integrator,
    platform,
    properties
  ) 
  simulation.context.setPositions(state.getPositions()) # This way we always start from the positions at the end of NPT equil.
  simulation.context.setVelocities(state.getVelocities())
  simulation.context.setPeriodicBoxVectors(
      *state.getPeriodicBoxVectors()
    )
  simulation.context.setVelocitiesToTemperature(md_temperature) # We generate different velocities using Boltzmann distribution to have different sim (we could add a seed to have the same, not needed).

  return simulation,system,integrator,pdb

def reporter_setup(simulation,
                   iter:int,
                   traj_name:str,
                   mdsteps:int,
                   reporter_interval:int        = 5_000,
                   traj_interval:int            = 1_000,
                   ):
  
  simulation.reporters.append(
    app.StateDataReporter(
      stdout,
      reporter_interval,
      totalSteps=mdsteps,
      step=True, potentialEnergy=True, temperature=True, volume=True, density=True, progress=True
    )
  )
  
  traj_filename = f'{traj_name}_{iter}.dcd'
  simulation.reporters.append(
    app.DCDReporter(traj_filename,traj_interval)
  )
  return simulation, traj_filename

def production_simulations( folder_prep:str              = 'system_prep',
                            timestep:unit.Quantity       = 2.0*unit.femtosecond,
                            friction_coeff:unit.Quantity = 1.0*1/unit.picosecond,
                            md_temperature:unit.Quantity = 298.15*unit.kelvin,
                            reporter_interval:int        = 5_000,
                            traj_interval:int            = 1_000,
                            traj_numb:int                = 3,
                            md_prod_total:int            = 150.0*unit.nanosecond,
                            ini:int                      = 2,
                            fin:int                      = 3
                            ):



  md_prod_length = md_prod_total / traj_numb
  len_md = str(md_prod_length.value_in_unit(unit.nanosecond))
  if len_md.split('.')[1] != '0':
    len_md = len_md.replace('.','_')
    dec_len_md = len_md.split('_')[1]
    print(len(dec_len_md))
    if len(dec_len_md) > 2:
      dec_len_md = dec_len_md[:2]
      len_md = len_md.split('_')[0] + '_' + dec_len_md
  else:
     len_md = len_md.split('.')[0]


  #* System and file names *#

  traj_name       = f'traj_{len_md}ns'

  for i in range(ini,fin+1):

    simulation,system,integrator,pdb = set_up_with_checkpoint(i,folder_prep)

    mdsteps = round(md_prod_length/timestep)

    simulation, traj_filename  = reporter_setup(simulation,
                                                i,
                                                traj_name,
                                                mdsteps,
                                                reporter_interval,
                                                traj_interval)

    print('\nRunning Production')

    simulation.step(mdsteps)

    print('\nDONE')

    save_checkpoint_for_restart(i,system,simulation,integrator)
  return traj_name