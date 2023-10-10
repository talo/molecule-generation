from io import StringIO, BytesIO
from typing import Any

from openbabel import pybel
from rdkit import Chem


RDKitMol = Any
PybelMol = Any


def rdkit_mol_to_sdf_string(rdkit_mol: RDKitMol) -> str:
    """
    Convert an RDKit mol to an sdf string.
    """
    sio = StringIO()
    writer = Chem.SDWriter(sio)
    writer.write(rdkit_mol)
    writer.close()
    sdf_string: str = sio.getvalue()
    return sdf_string


def pybel_mol_to_rdkit_mol(pybel_mol: PybelMol) -> RDKitMol:
    """
    Convert a pybel mol to an RDKit mol.
    """
    sdf_string: str = pybel_mol.write("sdf")
    sdf_bytes = BytesIO(sdf_string.encode("utf8"))
    rdkit_mol = next(Chem.ForwardSDMolSupplier(sdf_bytes, removeHs=False, sanitize=False))  # type: ignore
    return rdkit_mol


def pybel_add_Hs_to_rdkit_mol(rdkit_mol: RDKitMol) -> RDKitMol:
    """
    Add Hs to an RDKit mol using pybel.

    MUCH more robust than rdkit.Chem.AddH(mol, addCoords=True).
    """
    sdf_string = rdkit_mol_to_sdf_string(rdkit_mol)
    # convert to pybel mol
    pybel_mol = pybel.readstring("sdf", string=sdf_string)
    # add Hs with openbabel, this is MUCH more robust than Chem.AddH(mol, addCoords=True)
    pybel_mol.addh()
    # convert back to rdkit mol
    rdkit_mol_with_Hs = pybel_mol_to_rdkit_mol(pybel_mol)
    return rdkit_mol_with_Hs
