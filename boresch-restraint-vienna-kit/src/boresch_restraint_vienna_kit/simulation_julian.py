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
from sys import stdout
import numpy as np
import os, pathlib
from openmm.app import PDBFile,Modeller,ForceField
from pdbfixer import PDBFixer

def tester_for_bad_template(pdb_name:str):
  """ Tests the inputed PDB for OpenMM template errors.

  Args:
      pdb_name (str): name of the PDB, it puts the '_CLEAN.pdb' already from the workflow.
  """  
  print('Testing protein for bad template OpenMM error.')
  pdb = PDBFile(f'{pdb_name}_CLEAN.pdb')
  ff = ForceField('amber14/protein.ff14SB.xml')
  mod = Modeller(pdb.topology, pdb.positions)
  ff.createSystem(mod.topology)

def preparing_protein_from_PDBank(pdb_name:str,
                                  residues_to_remove:dict,
                                  protonate:bool):
  """Prepares the protein from a file obtained from the Protein Data Bank.

  Args:
      pdb_name (str): original pdb name from the beginning. It adds the '_fixed.pdb' inside.
      residues_to_remove (dict): residues that will be remove that are not from the protein, such as cofactors ions or ligands. No need to put the waters.
      protonate (bool): True if the protein needs to be protonated.
  """  
  print('Preparing Protein-only PDB:')
  pdb = PDBFile(f'{pdb_name}_fixed.pdb')

  modeller = Modeller(pdb.topology, pdb.positions)
  print(' - Deleting Waters')
  modeller.deleteWater()

  print(' - Deleting non-protein atoms.')
  atoms_to_delete = [
      atom for atom in modeller.topology.atoms()
      if atom.residue.name in residues_to_remove
  ]

  modeller.delete(atoms_to_delete)

  if protonate:
    print(' - Protonating protein.')
    forcefield = ForceField(
    	'amber/ff14SB.xml',
    	'amber/tip3p_standard.xml',
    	'amber/tip3p_HFE_multivalent.xml',
    	'amber/phosaa10.xml'
    )

    modeller.addHydrogens(forcefield=forcefield,
                          pH=7.4)

  print(' - Checking protein for connectivity problems.')
  for chain in modeller.topology.chains():
      for res1, res2 in zip(chain.residues(), list(chain.residues())[1:]):
          # Check if residues are sequential
          if res2.index != res1.index + 1:
              print('  - Chain break between:', res1, 'and', res2, '. NEEDS SOLVING!!!')

  print(f' - Saving clean protein-only PDB as: {pdb_name}_CLEAN.pdb')
  PDBFile.writeFile(modeller.topology, modeller.positions, open(f'{pdb_name}_CLEAN.pdb', 'w'))

  tester_for_bad_template(pdb_name)

def preparing_ligand_from_PDBank(sdf_orig_name:str,
                                lig_resname:str,
                                protonate:bool):
  """Prepares the ligand from a file obtained from the Protein Data Bank. It is important to check that the SDF downloaded from the Protein Data Bank is in the same coordinates than the ligand is in the PDB file. If not the protein will be elsewhere.

  Args:
      sdf_orig_name (str): sdf filename of the ligand without the extension '.sdf'
      lig_resname (str): resname of the ligand.
      protonate (bool): True if the ligand needs to be protonated. It does not check pH all atoms will be in their neutral states.
  """  
  supp = Chem.SDMolSupplier(f'{sdf_orig_name}.sdf', removeHs=False)
  mol = supp[0]

  if protonate:
    print(' - Protonating ligand.')
    mol = Chem.AddHs(mol, addCoords=True)

  print(f'Saving ligand-only SDF as lig_{lig_resname}.sdf')
  writer = Chem.SDWriter(f'lig_{lig_resname}.sdf')
  writer.write(mol)
  writer.close()
  sdf_path = f'lig_{lig_resname}.sdf'
  print('### If you dont want a NEUTRAL ligand change it manually in the SDF file!!!!\n')

  return

def preparing_from_PDBank(pdb_name:str,
                          lig_resname:str,
                          sdf_orig_name,
                          residues_to_remove:dict,
                          protonate:bool):
  """Prepares system obtained from the Protein Data Bank, we need the CIF of the system and the SDF of the ligand.

  Args:
      pdb_name (str): name of the CIF file without the extension '.cif'.
      lig_resname (str): resname of the ligand.
      sdf_orig_name (_type_): sdf filename of the ligand without the extension '.sdf'
      residues_to_remove (dict): residues that will be remove that are not from the protein, such as cofactors ions or ligands. No need to put the waters.
      protonate (bool): True if the system needs to be protonated. Usually the case if comes from XRAY spectroscopy.
  """  
  print('\nFixing PDB')
  fixer = PDBFixer(filename=f'{pdb_name}.cif')
  fixer.findMissingResidues()
  print(fixer.findMissingResidues)
  fixer.findMissingAtoms()
  print(fixer.findMissingAtoms)
  fixer.addMissingAtoms()

  print(f'Saving fixed PDB as {pdb_name}_fixed.pdb')
  with open(f'{pdb_name}_fixed.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

  preparing_protein_from_PDBank(pdb_name,
                                residues_to_remove,
                                protonate)
  
  preparing_ligand_from_PDBank(sdf_orig_name,
                                lig_resname,
                                protonate)
  return

def preparing_ligand(molecule_name:str):
  """Prepares the ligand for the simulation, assigning the partial charges to the atoms, to have a forcefield of the ligand.

  Args:
      molecule_name (str): name of the ligand sdf file, don't put the '.sdf' extension. 

  Returns:
      The parameters and the ligand to be used by OpenMM.
  """  
  #! Getting Ligand !#
  supp = Chem.SDMolSupplier(f'{molecule_name}.sdf', removeHs=False)
  orig_rdkit = supp[0]
  ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]

  #! Charging the Ligand !#
  # To have the ligand's FF parameters.

  charge_settings = OpenFFPartialChargeSettings(partial_charge_method='am1bcc', off_toolkit_backend='ambertools')

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
                        forcefield:list = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'],
                        forcefield_small_molecule:str = 'openff-2.2.1' 
                        ):
  """Does the preparation of a system that has already been solvated. It assigns the cell dimensions for the PBC.

  Args:
      pdb_name (str): filename of the PDB with the solvated system.
      orig_rdkit (_type_): information of the ligand produced by the preparing function above.
      forcefield (list, optional): list of forcefields to be used. Defaults to ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'].
      forcefield_small_molecule (str, optional): forcefield to be used for the ligand. Defaults to 'openff-2.2.1'.

  Returns:
      system, topology, positions to be used for the simulation.
  """    
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

def solvating_system(pdb_name:str,
                     ligands,
                     forcefield:list = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'],
                     forcefield_small_molecule:str = 'openff-2.2.1',
                     protein_ligand_pdb = 'protein_ligand_system.pdb',
                     solvated_pdb = 'solvated_system.pdb',
):
  """Solvates a system where we only have the protein and the ligand.

  Args:
      pdb_name (str): filename of the PDB with the protein.
      ligands : we get this parameter from the preparing_ligand function.
      forcefield (list, optional): list of forcefields to be used. Defaults to ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'].
      forcefield_small_molecule (str, optional): forcefield to be used for the ligand. Defaults to 'openff-2.2.1'.
      protein_ligand_pdb (str, optional): name of the PDB to be written of the ligand and protein. Defaults to 'protein_ligand_system.pdb'.
      solvated_pdb (str, optional): name of the PDB to be written of the solvated system. Defaults to 'solvated_system.pdb'.

  Returns:
      system, topology, positions for the OpenMM simulation
  """  
  protein = app.PDBFile(pdb_name)

  off_mol = ligands[0].to_openff()
  ligand_top = off_mol.to_topology().to_openmm()
  ligand_pos = off_mol.conformers[0].to_openmm()

  modeller = app.Modeller(protein.topology, protein.positions)
  modeller.add(ligand_top,ligand_pos)

  with open(protein_ligand_pdb, 'w') as f:
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
    padding=1.5*unit.nanometer,
    ionicStrength = 0.15 * unit.molar
  )

  with open(solvated_pdb, 'w') as f:
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
  """We generate the integrator and simulation with the given parameters.

  Args:
      md_temperature (unit.Quantity): temperature of the simulation.
      friction_coeff (unit.Quantity): friction coefficient for the simulation.
      timestep (unit.Quantity): timestep in femtoseconds.
      topology : topology of the system.
      system : system variable with properties
      positions : positions of the atoms in the system 

  Returns:
      Integrator and simulation variables from OpenMM.
  """  
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
                 min_tolerance: unit.Quantity,
                 minimization_max_steps:int,
                 pdb_minimization_final:str
                 ):
  """This functions performs a minimization of the given simulation.

  Args:
      simulation : simulation, can be obtained with function "setting_up_simulation"
      min_tolerance (unit.Quantity): minimum energy tolerance per unit of length.
      minimization_max_steps (int): maximum number of minimization steps.
      pdb_minimization_final (str): namefile for the PDB generated after the minimization.

  Returns:
      simulation after minimization.
  """  
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
                   md_equil_nvt_length:unit.Quantity = 0.25*unit.nanosecond,
                   timestep:unit.Quantity       = 2.0*unit.femtosecond,                   
                   md_temperature:unit.Quantity = 298.15*unit.kelvin,
                   temp_increment:int = 100,
                   md_steps_after_ramp:int = 0,
                   pdb_heating_final:str = 'NVT_system.pdb'
                   ):
  """Performs a heating NVT simulation, it performs the amount of stages of heating for a given temperature increment in a given length of simulation.
  Args:
      simulation: simulation, can be obtained with function "setting_up_simulation"
      integrator: integrator, can be obtained with function "setting_up_simulation"
      reporter_interval (int): every this number of steps the reporter outputs information. 
      md_equil_nvt_length (unit.Quantity, optional): length of the NVT simulation in units of times such as nanoseconds.. Defaults to 0.25*unit.nanosecond.
      timestep (unit.Quantity, optional): timesteps of the simulation, must be the same that was given to the integrator. Defaults to 2.0*unit.femtosecond.
      md_temperature (unit.Quantity, optional): Temperature of the simulation, objective temperature for the ramp.. Defaults to 298.15*unit.kelvin.
      temp_increment (int, optional): Increments of temperature, the NVT simulation will be divided into stages at different temperatures following this increment. Defaults to 100.
      md_steps_after_ramp (int, optional): number of steps that want to be done after finishing with the NVT simulation, to let the system stabilize further. Defaults to 0.
      pdb_heating_final (str, optional): Name of the PDB file to output the system after the NVT simulation.  Defaults to 'NVT_system.pdb'. Defaults to 'NVT_system.pdb'.

  Returns:
      simulation after NVT equilibration
  """  
 
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
                   md_temperature:unit.Quantity = 298.15*unit.kelvin,
                   md_equil_npt_length:unit.Quantity = 0.5*unit.nanosecond,
                   md_pressure:unit.Quantity = 1*unit.bar,
                   pdb_density_equil_final:str = 'NPT_system.pdb'
                   ):
  """We perform an NPT simulation, at the beginning we add the barostat.

  Args:
      simulation : simulation, can be obtained with function "setting_up_simulation"
      system : OpenMM system were the barostat is added. 
      timestep (unit.Quantity, optional): timestep in units of times, usually femtosecond. Defaults to 2.0*unit.femtosecond.
      md_temperature (unit.Quantity, optional): Temperature of the simulation, usually in Kelvin. Defaults to 298.15*unit.kelvin.
      md_equil_npt_length (unit.Quantity, optional): Length of the NPT simulation, usually in nanoseconds. Defaults to 0.5*unit.nanosecond.
      md_pressure (unit.Quantity, optional): Pressure of the simulation, usually in bar. Defaults to 1*unit.bar.
      pdb_density_equil_final (str, optional): Name of the PDB file to output the system after the NVT simulation. Defaults to 'NPT_system.pdb'.

  Returns:
      simulation after NPT equilibration, system with barostat, npt_state to be able to start several simulations from the same last step of the NPT equilibration.
  """  

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

def save_checkpoint_for_restart(iter:int,system,simulation,integrator):
  """We save several files to use as checkpoints, to be able to restart simulations from them.

  Args:
      iter (int): which step of the simulation is this, I usually use 0 for preparation simulation and then add 1 for each simulation we have done.
      system : OpenMM system.
      simulation: simulation, can be obtained with function "setting_up_simulation"
      integrator: integrator, can be obtained with function "setting_up_simulation"

  Returns:
      Several xml with information on the system, state and integrator, useful when we run in a different machine; and a checkpoint file, from which a system can be restarted, we only need this one when we are on the same machine.
  """  
  with open(f'system_{iter}.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))
  with open(f'integrator_{iter}.xml', 'w') as f:
    f.write(XmlSerializer.serialize(integrator))
  state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=False,
    enforcePeriodicBox=True
  )
  with open(f'state_{iter}.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))

  simulation.saveCheckpoint(f'checkpoint_{iter}.chk')
  positions = state.getPositions()

  with open(f'last_frame_{iter}.pdb', 'w') as f:
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
  """It does all the preparation steps but the preparing from PDBank functions. It generates the system, minimizes, NVT equilibration and NPT equilibration.

  Args:
      molecule_name (str): name of the ligand sdf file, don't put the '.sdf' extension.
      pdb_name (str): filename of the PDB with the solvated system (solvated_PDB = True) or filename of the PDB with the protein (solvated_PDB = False). Put the extension '.pdb' 
      solvated_PDB (bool): True if pdb name is already solvated.
      folder_prep (str, optional): name of the folder where all the preparatory steps will be done. Defaults to 'system_prep'.
      timestep (unit.Quantity, optional): timestep in units of times, usually femtosecond. Defaults to 2.0*unit.femtosecond.
      friction_coeff (unit.Quantity, optional): friction_coeff (unit.Quantity): friction coefficient for the simulation.
      md_temperature (unit.Quantity, optional): temperature at which you want the simulation, usually in kelvin. Defaults to 298.15*unit.kelvin.
      minimization_max_steps (int, optional): maximum number of minimization steps. Defaults to 5000.
      min_tolerance (unit.Quantity, optional): minimum energy tolerance per unit of length. Defaults to 4*unit.kilojoule_per_mole/unit.nanometer.
      pdb_minimization_final (str, optional): namefile for the PDB generated after the minimization. Defaults to 'minimized_system.pdb'.
      reporter_interval (int, optional): every this number of steps the reporter outputs information. Defaults to 5_000.
      temp_increment (int, optional): Increments of temperature, the NVT simulation will be divided into stages at different temperatures following this increment. Defaults to 100.
      md_steps_after_ramp (unit.Quantity, optional): number of steps that want to be done after finishing with the NVT simulation, to let the system stabilize further. Defaults to 0.
      md_equil_npt_length (unit.Quantity, optional): Length of the NPT simulation, usually in nanoseconds. Defaults to 0.5*unit.nanosecond.
      md_pressure (unit.Quantity, optional): Pressure of the simulation, usually in bar. Defaults to 1*unit.bar.
      pdb_NPT_final (str, optional): Name of the PDB file to output the system after the NVT simulation. Defaults to 'NPT_system.pdb'.
  """  

  path = pathlib.Path(f'./{folder_prep}')
  path.mkdir(exist_ok=True)
  os.chdir(folder_prep)

  os.system(f'cp ../{pdb_name} .')
  os.system(f'cp ../{molecule_name}.sdf .')

  orig_rdkit,ligands = preparing_ligand(molecule_name)
  
  if solvated_PDB:
    system, topology, positions = system_solvated_PDB(pdb_name,
                                                    orig_rdkit,
                                                    forcefield = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'],
                                                    forcefield_small_molecule = 'openff-2.2.1' )

  else:
    system, topology, positions = solvating_system(pdb_name,
                                                 ligands,
                                                 forcefield = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml', 'amber/phosaa10.xml'],
                                                 forcefield_small_molecule = 'openff-2.2.1',
                                                 protein_ligand_pdb = 'protein_ligand_system.pdb',
                                                 solvated_pdb = 'solvated_system.pdb',)

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
                                               md_temperature,
                                               md_equil_npt_length,
                                               md_pressure,
                                               pdb_NPT_final)
  
  save_checkpoint_for_restart(0,system,simulation,integrator)

  os.chdir('..')

def set_up_with_checkpoint(iter:int,
                           folder_prep:str,
                           timestep:unit.Quantity       = 2.0*unit.femtosecond,
                           friction_coeff:unit.Quantity = 1.0*1/unit.picosecond,
                           md_temperature:unit.Quantity = 298.15*unit.kelvin,
):
  """Sets up a simulation to restart one finished previously.

  Args:
      iter (int): index of the simulation about to start, not the one with the checkpoint.
      folder_prep (str): foldername of the directory with the preparation simulations.
      timestep (unit.Quantity, optional): timestep in units of times, usually femtosecond. Defaults to 2.0*unit.femtosecond.
      friction_coeff (unit.Quantity, optional): friction coefficient for the simulation. Defaults to 1.0*1/unit.picosecond.
      md_temperature (unit.Quantity, optional): temperature of the simulation. Defaults to 298.15*unit.kelvin.

  Returns:
      simulation, system, integrator and pdb information from the previous simulation.
  """  
  
  if iter == 1:
    folder = folder_prep
  else:
    folder = '.'
  
  pdb = app.PDBFile(f'{folder}/last_frame_{iter-1}.pdb')

  with open(f'{folder}/system_{iter-1}.xml') as f:
      system = XmlSerializer.deserialize(f.read())

  with open(f'{folder}/integrator_{iter-1}.xml') as f:
      integrator = XmlSerializer.deserialize(f.read())

  platform = Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'mixed'}  # MUST match original run
  simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

  with open(f'{folder}/checkpoint_{iter-1}.chk', 'rb') as f:
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
  """Generates the reporters for the output and the trajectory.

  Args:
      simulation : simulation, usually obtained from 'set_up_with_checkpoint' function. 
      iter (int): index of the simulation to do.
      traj_name (str): trajectory name before the index and the '.dcd' extension.
      mdsteps (int): total amount of steps to do, usually given previously as unit of time and converted automatically with time_step.
      reporter_interval (int, optional): the amount of steps before reporting again in the output. Defaults to 5_000.
      traj_interval (int, optional): the amount of steps before writing trajectory step again. Defaults to 1_000.

  Returns:
      simulation with the reporters, traj_filename
  """  
  
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
                            ini:int                      = 1,
                            fin:int                      = 3
                            ):
  """Runs the production simulations restarting from the preparation simulation (ini = 1) or from another production trajectory (i > 1).

  Args:
      folder_prep (str, optional): foldername of the directory with the preparation simulations. Defaults to 'system_prep'.
      timestep (unit.Quantity, optional): timestep in units of times, usually femtosecond. Defaults to 2.0*unit.femtosecond.
      friction_coeff (unit.Quantity, optional): friction coefficient for the simulation. Defaults to 1.0*1/unit.picosecond.
      md_temperature (unit.Quantity, optional): temperature of the simulation. Defaults to 298.15*unit.kelvin.
      reporter_interval (int, optional): the amount of steps before reporting again in the output. Defaults to 5_000.
      traj_interval (int, optional): the amount of steps before writing trajectory step again. Defaults to 1_000.
      traj_numb (int, optional): Amount of trajectories you want to divide the md_prod_total in. Defaults to 3.
      md_prod_total (int, optional): Total lenght in units of time of the simulation, usually in nanosecond. Defaults to 150.0*unit.nanosecond.
      ini (int, optional): index of the first trajectory to do. Defaults to 1.
      fin (int, optional): index of the last trajectory to do. Defaults to 3.

  Returns:
      traj_name, filename before the index and the '.dcd' extension.
  """  

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

    simulation,system,integrator,pdb = set_up_with_checkpoint(i,folder_prep,timestep,friction_coeff,md_temperature)

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