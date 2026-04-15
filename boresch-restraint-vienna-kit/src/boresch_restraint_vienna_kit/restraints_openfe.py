from openmm import *

from openfe.protocols.openmm_afe.equil_afe_settings import BoreschRestraintSettings
from openfe.protocols.openmm_afe.abfe_units import ABFEComplexSetupUnit
from openff.units import unit as off_unit
from openmm.app import PDBFile
from pathlib import Path
from rdkit import Chem
from pint import UnitRegistry


def restraint_search_openfe(pdb_name_search: str,
                            guest_sdf_name: str,
                            guest_resname: str,
                            trajectory_name: str,
                            path :str = '.',
                            temperature :float = 298.15 
):

  base = Path(path)
  pdb = PDBFile(pdb_name_search)
  supp = Chem.SDMolSupplier(guest_sdf_name, removeHs=False)

  guest_rdmol = supp[0]

  universe_mda = ABFEComplexSetupUnit._get_mda_universe(
    topology = pdb.topology,
    positions =  None,
    trajectory = base / trajectory_name,
    )

  guest = universe_mda.select_atoms(f"resname {guest_resname}")
  host = universe_mda.select_atoms("protein")

  guest_atom_ids = guest.indices.astype(int).tolist()
  host_atom_ids = host.indices.astype(int).tolist()


  settings = BoreschRestraintSettings()
  temperature = temperature * off_unit('kelvin')

  geom, restrain = ABFEComplexSetupUnit._get_boresch_restraint(
    universe = universe_mda,
    guest_rdmol = guest_rdmol,
    guest_atom_ids = guest_atom_ids,
    host_atom_ids = host_atom_ids,
    temperature = temperature,
    settings = settings,
  )

  ureg = UnitRegistry()

  print('\n=== "OPENFE RESULTS" ===')
  print('\nGUEST ATOMS')
  print(geom.guest_atoms)
  print('HOST ATOMS')
  print(geom.host_atoms)
  print('\n EQUILIBRIUM PARAMETERS')
  print("r_aA0:", geom.r_aA0.to(ureg.angstrom))
  print("theta_A0:", geom.theta_A0.to(ureg.degree))
  print("theta_B0:", geom.theta_B0.to(ureg.degree))
  print("phi_A0:", geom.phi_A0.to(ureg.degree))
  print("phi_B0:", geom.phi_B0.to(ureg.degree))
  print("phi_C0:", geom.phi_C0.to(ureg.degree))
  # print('\nALL')
  # print(geom)

  last_frame_vars = [float(geom.r_aA0.to("angstrom").magnitude), float(geom.theta_A0.to("degree").magnitude), float(geom.theta_B0.to("degree").magnitude), float(geom.phi_A0.to("degree").magnitude), float(geom.phi_B0.to("degree").magnitude), float(geom.phi_C0.to("degree").magnitude)]

  return universe_mda, geom.host_atoms, geom.guest_atoms, last_frame_vars