from rdkit import Chem
from rdkit.Chem import Draw,rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def williams_anchors_to_names(universe_mda, anchors):
  """We convert William's anchors results to obtain the resnames and resids of LIG and anchored residue, also we will convert the index of the anchors into their atom names.

  Args:
      universe_mda (MDAnalysis universe): MDAnalysis universe from the William's procedure.
      anchors (list): atom indexes of the different anchors found by William's procedure.

  Returns:
      resid and resname of the Ligand,resid and resname of the anchored residue, and a list of the atom names of the anchors.
  """  
  h0 = universe_mda.atoms[anchors[0]]
  h1 = universe_mda.atoms[anchors[1]]
  h2 = universe_mda.atoms[anchors[2]]
  g0 = universe_mda.atoms[anchors[3]]
  g1 = universe_mda.atoms[anchors[4]]
  g2 = universe_mda.atoms[anchors[5]]

  resname_L = g0.resname
  resname_P = h0.resname

  resid_L = g0.resid
  resid_P = h0.resid

  atom_names = [h0.name,h1.name,h2.name,g0.name,g1.name,g2.name]

  return resname_L, resname_P, resid_L, resid_P, atom_names

def openfe_anchors_to_names(universe_mda, host_atoms,guest_atoms):
  """We convert OpenFE's anchors results to obtain the resnames and resids of LIG and anchored residue, also we will convert the index of the anchors into their atom names. It is important to use the MDA universe of William's procedure, because OpenFE changes the resid and makes it different from the PDB.

  Args:
      universe_mda (MDAnalysis universe): MDAnalysis universe from the William's procedure.
      host_atoms (list): atom indexes of the different anchors of the protein found by OpenFE's procedure.
      guest_atoms (list): atom indexes of the different anchors of the ligand found by OpenFE's procedure.

  Returns:
      resid and resname of the Ligand,resid and resname of the anchored residue, and a list of the atom names of the anchors.
  """  ""
  h0 = universe_mda.atoms[host_atoms[0]]
  h1 = universe_mda.atoms[host_atoms[1]]
  h2 = universe_mda.atoms[host_atoms[2]]
  g0 = universe_mda.atoms[guest_atoms[0]]
  g1 = universe_mda.atoms[guest_atoms[1]]
  g2 = universe_mda.atoms[guest_atoms[2]]

  resname_L = g0.resname
  resname_P = h0.resname

  resid_L = g0.resid
  resid_P = h0.resid

  atom_names = [h0.name,h1.name,h2.name,g0.name,g1.name,g2.name]

  return resname_L, resname_P, resid_L, resid_P, atom_names

def load_pdb_without_waters(path_pdb):
    with open(path_pdb, "r") as f:
        lines = f.readlines()

    filtered_lines = [
        line for line in lines
        if not (line.startswith("HETATM") and "HOH" in line[17:21])
        and not (line.startswith("HETATM") and "WAT" in line[17:21])
        and not (line.startswith("HETATM") and "SOL" in line[17:21])
    ]

    pdb_block = "".join(filtered_lines)

    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("RDKit failed to parse filtered PDB")

    return mol

def extract_submol(mol, selection):
  """Return a new Mol containing only selected atoms"""
  emol = Chem.RWMol()

  idx_map = {}

  for atom in mol.GetAtoms():
    info = atom.GetPDBResidueInfo()
    if not info:
      continue

    key = (info.GetResidueName().strip(), info.GetResidueNumber())

    if key in selection:
      new_idx = emol.AddAtom(atom)
      idx_map[atom.GetIdx()] = new_idx

  # add bonds
  for bond in mol.GetBonds():
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()

    if a1 in idx_map and a2 in idx_map:
      emol.AddBond(
          idx_map[a1],
          idx_map[a2],
          bond.GetBondType()
      )

  return emol.GetMol(), idx_map

def find_atom(mol, atomname):
  for atom in mol.GetAtoms():
    info = atom.GetPDBResidueInfo()
    if info and info.GetName().strip() == atomname:
      return atom.GetIdx()
  return None

def drawing(path_pdb:str, 
            path_lig_sdf:str, 
            resname_lig:str,
            resid_lig :int, 
            resname_prot:str,
            resid_prot:int,
            anchors_names:list,
            figure_name:str
            ):
  # mol = Chem.MolFromPDBFile(path_pdb, removeHs=False)
  mol = load_pdb_without_waters(path_pdb)

  lig_sdf = Chem.SDMolSupplier(path_lig_sdf)[0]

  lig_pdb, lig_map = extract_submol(mol, {(resname_lig, resid_lig)})

  pdb_atoms = list(lig_pdb.GetAtoms())
  sdf_atoms = list(lig_sdf.GetAtoms())

  atom_map = {
      pdb_atoms[i].GetIdx(): sdf_atoms[i].GetIdx()
      for i in range(min(len(pdb_atoms), len(sdf_atoms)))
  }

  lig = lig_sdf

  res, res_map = extract_submol(mol, {(resname_prot, resid_prot)})

  h0 = find_atom(res, anchors_names[0])
  h1 = find_atom(res, anchors_names[1])
  h2 = find_atom(res, anchors_names[2])

  g0_pdb = find_atom(lig_pdb, anchors_names[3])
  g0 = atom_map[g0_pdb]
  g1_pdb = find_atom(lig_pdb, anchors_names[4])
  g1 = atom_map[g1_pdb]
  g2_pdb = find_atom(lig_pdb, anchors_names[5])
  g2 = atom_map[g2_pdb]

  combo = Chem.CombineMols(lig, res)
  ed_combo = Chem.EditableMol(combo)

  offset = lig.GetNumAtoms()

  if anchors_names[0] == 'N':

    C = find_atom(res, 'C')
    anchor_C = C + offset

    dummy = ed_combo.AddAtom(Chem.Atom(0))

    bond_type = Chem.BondType.HYDROGEN

    ed_combo.AddBond(anchor_C, dummy, bond_type)
  else: 
    N = find_atom(res, 'N')
    anchor_N = N + offset

    dummy = ed_combo.AddAtom(Chem.Atom(0))

    bond_type = Chem.BondType.HYDROGEN

    ed_combo.AddBond(anchor_N, dummy, bond_type)


  ed_combo.AddBond(
      g0,
      h0 + offset,
      bond_type
  )

  newmol = ed_combo.GetMol()

  newmol.GetAtomWithIdx(dummy).SetProp(
    "atomNote",
    f"{resname_prot}{resid_prot}"
  )

  highlight_atoms = [g0,g1,g2, h0 + offset,h1 + offset,h2 + offset]#, dummy]
  highlight_colors = {
      g0: (0.70, 0.96, 0.79),       
      g1: (0.99, 0.92, 0.55),        
      g2: (0.62, 0.81, 1.00),       
      h0 + offset: (0.70, 0.96, 0.79), 
      h1 + offset: (0.99, 0.92, 0.55), 
      h2 + offset: (0.62, 0.81, 1.00), 
      # dummy : (0.01, 0.5, 1)
  }


  rdDepictor.Compute2DCoords(newmol)

  drawer = rdMolDraw2D.MolDraw2DSVG(1200,800)
  drawer.drawOptions().explicitMethyl = True
  drawer.drawOptions().atomHighlightsAreCircles = True
  drawer.drawOptions().bondLineWidth = 5.0
  drawer.drawOptions().annotationFontScale = 1.1
  # drawer.drawOptions().annotationColour = (0.01, 0.5, 1)

  drawer.DrawMolecule(
                      newmol,
                      highlightAtoms=list(highlight_atoms),
                      highlightBonds=[],
                      highlightAtomColors=highlight_colors,
  )
  drawer.FinishDrawing()

  svg = drawer.GetDrawingText()

  with open(f'{figure_name}.svg', 'w') as f:
    f.write(svg)