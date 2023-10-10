import os
import random
import string
import subprocess
import time
from functools import wraps, partial
from io import StringIO, BytesIO
from pathlib import Path
from typing import Any, Callable

import numpy as np
import numpy.typing as npt
import pandas as pd
from loguru import logger
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

from molecule_generation import load_model_from_directory
from mso.optimizer import BasePSOptimizer
from mso.objectives.scoring import ScoringFunction


RDKitMol = Any
PybelMol = Any
Score = float


def check_valid_mol(func: Callable[[Any], Any]) -> Callable[[Any], Any]:
    """
    Decorator function that checks if a mol object is None (resulting from a non-processable SMILES string)
    :param func: the function to decorate.
    :return: The decorated function.
    """
    @wraps(func)
    def wrapper(mol, *args, **kwargs):
        if mol is not None:
            return func(mol, *args, **kwargs)
        else:
            return 0
    return wrapper
    

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
    # NOTE: do need to do sanitize for 3d embed
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


def run_gnina_single(
    mol_with_Hs: RDKitMol,
    protein_path: Path,
    autobox_ligand_path: Path,
    lig_name: str = "gnina",
    exhaustiveness: int = 64,
    autobox_add: int = 10,
    num_poses: int = 1,  # or 10
    num_cpus: int = 24,
    seed: int = 1337,
    temp_dir: Path = Path("./moler_oracle_temp/")
) -> Path:
    """
    Run gnina scoring over a single molecule.

    NOTE: needs about 3 GB of GPU VRAM.
    """
    temp_dir.mkdir(exist_ok=True, parents=True)

    # Write mol with Hs to .sdf
    random_hash = ''.join(
        random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(20)
    )
    output_mol_with_Hs_path = temp_dir / f"docked_in_{lig_name}_{random_hash}.sdf"
    output_gnina_path = temp_dir / f"docked_out_{lig_name}_{random_hash}.sdf"

    if output_mol_with_Hs_path.exists():
        print(f"output_mol_with_Hs_path {output_mol_with_Hs_path} already exists, deleting it.")
        os.remove(output_mol_with_Hs_path)
    if output_gnina_path.exists():
        print(f"output_gnina_path {output_gnina_path} already exists, deleting it.")
        os.remove(output_gnina_path)

    writer = Chem.SDWriter(str(output_mol_with_Hs_path))
    writer.write(mol_with_Hs)
    writer.close()

    command = (
        f"gnina --receptor {protein_path} "
        f"--ligand {output_mol_with_Hs_path} "
        f"--out {output_gnina_path} "
        f"--autobox_ligand {autobox_ligand_path} "
        f"--exhaustiveness {exhaustiveness} "
        f"--autobox_add {autobox_add} "
        f"--num_modes {num_poses} "  # we usually only need 1 pose
        f"--seed {seed} "
        f"--cpu {num_cpus}"
    )
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    process.communicate()
    return output_gnina_path


def combine_ligand_with_protein_pdb(protein_pdb_path: Path, lig_rdkit_mol: RDKitMol, output_path: Path):
    # make ligand PDB from rdkit mol
    top_pose_with_Hs = pybel_add_Hs_to_rdkit_mol(lig_rdkit_mol)
    Chem.SanitizeMol(top_pose_with_Hs)
    
    # NOTE: it is safer to go via SDF, preserving all ligand connectivity info
    sdf_string = rdkit_mol_to_sdf_string(top_pose_with_Hs)
    pybel_mol = pybel.readstring("sdf", string=sdf_string)
    lig_pdb_string: str = pybel_mol.write("pdb")
    
    # get lig PDB lines
    lig_pdb_lines = [
        line for line in lig_pdb_string.split("\n") 
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT")
    ]

    # load protein PDB lines
    with open(protein_pdb_path, "r") as f:
        prot_pdb_lines = [line.strip() for line in f.readlines()]
        prot_pdb_lines = [line for line in prot_pdb_lines if line.startswith("ATOM")]

    # combine lines
    prot_lig_pdb_lines = prot_pdb_lines.copy()
    prot_lig_pdb_lines.extend(lig_pdb_lines)

    # save combined PDB with openbabel's pybel
    pybel_prot_lig_pdb = pybel.readstring("pdb", string="\n".join(prot_lig_pdb_lines))
    pybel_prot_lig_pdb.write("pdb", str(output_path), overwrite=True)


def retrieve_plip_interactions(complex_pdb_path: Path):
    """
    Retrieves the interactions from PLIP.

    Parameters
    ----------
    complex_pdb_path :
        Path to PDB file of the complex.

    Returns
    -------
    dict :
        A dictionary mapping all possible interactions types to interacting residues.
    """
    complex = PDBComplex()
    complex.load_pdb(str(complex_pdb_path))  # load the pdb file

    assert len(complex.ligands) == 1, "Error, only 1 ligand is allowed!"
    for ligand in complex.ligands:
        complex.characterize_complex(ligand)  # find ligands and analyze interactions
        
    # for key, site in sorted(complex.interaction_sets.items()):
    site_key = list(complex.interaction_sets.keys())[0]
    site = complex.interaction_sets[site_key]
    binding_site = BindingSiteReport(site)  # collect data about interactions
    # tuples of *_features and *_info will be converted to pandas DataFrame
    keys = (
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    )
    # interactions is a dictionary which contains relevant information for each
    # of the possible interactions: hydrophobic, hbond, etc. in the considered
    # binding site. Each interaction contains a list with
    # 1. the features of that interaction, e.g. for hydrophobic:
    # ('RESNR', 'RESTYPE', ..., 'LIGCOO', 'PROTCOO')
    # 2. information for each of these features, e.g. for hydrophobic
    # (residue nb, residue type,..., ligand atom 3D coord., protein atom 3D coord.)
    return {
        k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")
        for k in keys
    }
    

@check_valid_mol
def dock_and_get_interactions(
    mol: RDKitMol, 
    protein_pdb_path: Path,
    autobox_ligand_path: Path,
    residues_of_interest: set[str],
    necessary_residues: set[str],
    exhaustiveness: int = 64,
    seed: int = 1337,
    num_poses: int = 10,
    temp_dir: Path = Path("./moler_oracle_temp/"),
    debug: bool = False
) -> Score:
    assert len(residues_of_interest) > 0
    assert len(necessary_residues) > 0
    assert len(necessary_residues & residues_of_interest) == len(necessary_residues)
    
    # always add Hs
    mol_with_Hs = pybel_add_Hs_to_rdkit_mol(mol)
    
    # take mol, embed RDKit conformer
    # if RDKit fails to get a conformer, then return a score of 0
    mol_with_Hs.UpdatePropertyCache()  # needed before embed
    embed_status = AllChem.EmbedMolecule(mol_with_Hs, randomSeed=seed, maxAttempts=500)
    if embed_status == -1:  # failed, try again with random coords
        params = AllChem.ETKDGv2()
        params.useRandomCoords = True
        embed_status = AllChem.EmbedMolecule(mol_with_Hs, params, randomSeed=seed, maxAttempts=500)
        if embed_status == -1:  # failed again, return score of 0
            return 0
    
    # then dock with GNINA
    if debug: logger.debug("began docking")
    docked_poses_path = run_gnina_single(
        mol_with_Hs,
        protein_pdb_path,
        autobox_ligand_path=autobox_ligand_path,
        exhaustiveness=exhaustiveness,
        num_poses=num_poses,
        seed=seed,
        temp_dir=temp_dir
    )
    hash = docked_poses_path.stem.replace("docked_out_gnina_", "")
    if debug: logger.debug("finished docking")
    
    # load docked sdf + protein PDB with plip, get interactions
    # N poses --> get N complexes --> analyse all by PLIP
    # checking all N poses helps eliminate variance from randomness of RDKit conformer generator & GNINA's finite exhaustiveness
    docked_poses = Chem.SDMolSupplier(str(docked_poses_path))

    interaction_score_per_pose = []
    best_pose = None
    best_pose_score = float("-inf")
    best_pose_idx = -1
    best_residues = set()
    for i, pose in enumerate(docked_poses):
        prot_lig_pdb_path = docked_poses_path.parent / f"{docked_poses_path.stem}_complex_pose{i}.pdb"
        combine_ligand_with_protein_pdb(protein_pdb_path, pose, prot_lig_pdb_path)

        # analyse with PLIP
        # https://projects.volkamerlab.org/teachopencadd/talktorials/T016_protein_ligand_interactions.html
        # Profiling-protein-ligand-interactions-using-PLIP
        inters_to_residues = retrieve_plip_interactions(prot_lig_pdb_path)
        all_residues = set()
        for inter_type, residues in inters_to_residues.items():
            residues_set = {f"{res[1]}{res[0]}" for res in residues[1:]}
            all_residues |= residues_set

        if len(necessary_residues & all_residues) < len(necessary_residues):
            pose_score = 0
        else:
            pose_score = len(residues_of_interest & all_residues) / len(residues_of_interest)
        if pose_score > best_pose_score:
            best_pose_score = pose_score
            best_pose = pose
            best_pose_idx = i
            best_residues = all_residues
        interaction_score_per_pose.append(pose_score)
        
        # cleanup
        os.remove(prot_lig_pdb_path)

    # add list of interacting residues as property & update docked out SDF --> resave with float score in name too
    best_pose.SetProp("residues", ",".join(sorted(list(best_residues), key=lambda x: int(x[3:]))))
    best_pose.SetProp("smiles", Chem.MolToSmiles(best_pose))
    best_pose_out_path = docked_poses_path.parent / f"best_score{best_pose_score:.4f}_pose{best_pose_idx}_{hash}.sdf"
    sdf_writer = Chem.SDWriter(str(best_pose_out_path))
    sdf_writer.write(best_pose)
    
    return best_pose_score


class MoLeRModel:
    def __init__(self, model, scaffold: str, debug: bool = False):
        self.model = model
        self.scaffold = scaffold
        self.debug = debug
        
    def emb_to_seq(self, x: npt.NDArray[np.float32]) -> list[str]:
        # x is 2D np array of shape [num_particles, emb_dim]
        if self.debug: logger.debug(f"shape of x in emb_to_seq: {x.shape}")
        decoded_smiles = self.model.decode(x, scaffolds=[self.scaffold] * x.shape[0])
        return decoded_smiles

    def seq_to_emb(self, smiles: list[str]) -> npt.NDArray[np.float32]:
        if self.debug: logger.debug(f"smiles in seq_to_emb: {smiles}")
        if isinstance(smiles, str):
            embeddings = np.array(self.model.encode([smiles]))
        else:
            embeddings = np.array(self.model.encode(smiles))
        if self.debug: logger.debug(f"embeddings in seq_to_emb: {embeddings.shape}")
        return embeddings


def ic50_str_to_float(ic50_str: str) -> float:
    if isinstance(ic50_str, float):
        return ic50_str
    elif ">" in ic50_str:
        return float(ic50_str.replace(">", ""))
    else:
        return float(ic50_str)
    

def optimise():
    np.random.seed(1337)
    random.seed(1337)
    os.environ["PYTHONHASHSEED"] = str(1337)

    time_hash = int(time.time())

    # run_name = "50particles_5rounds_test_v1"  # PRT1009303
    # run_name = "100parts_10steps_5res_5poses_8exh_v1_PRT1008447"
    run_name = "TEST_100parts_10steps_5res_5poses_8exh_v1_top100PRTs"
    num_parts = 3

    # smiles & scaffold to initialise
    scaffold = "CNC(=O)C1=C(NC2=CC=CC=C2OC)C=C(NC(=O)C2CC2)N=N1"

    xl_file = pd.read_excel("/home/minhtoo/molecule-generation/data/JAK2_IC50_data_Sep18.xlsx", sheet_name=None)
    df = xl_file["browser export"]
    df["IC50_NM_JH2_V617F"] = df["IC50_NM_JH2_V617F"].map(ic50_str_to_float)
    df = df.sort_values(by="IC50_NM_JH2_V617F")  # ascending IC50, from lowest (strongest binder) to highest (weakest binder)
    ref_smis = df["SMILES"].tolist()
    ref_smis = [smi.replace("|r|", "").strip() for smi in ref_smis]  # deal with artifacts in excel
    ref_smis = [smi for smi in ref_smis if len(smi) < 120]  # some SMILES are very long, we prolly dont want them for now
    init_smiles = ref_smis[:num_parts]
    # init_smiles = "CNC(=O)c1nnc(cc1Nc1cccc(c1OC)-c1ccc(nc1)C(=O)N(C)Cc1cc(ccn1)-c1ccncc1)NC(=O)C1CC1"  # PRT1009303
    # init_smiles = "CNC(=O)c1nnc(cc1Nc1cccc(c1OC)-c1cnn(c1)C1CN(C1)Cc1ccccn1)NC(=O)C1CC1"  # PRT1008447

    # define oracle
    protein_pdb_path = Path("/home/minhtoo/molecule-generation/data/mutant_PRT1008596.align.fixed.addH.pdb")
    autobox_ligand_path = Path("/home/minhtoo/molecule-generation/data/PRT1008596.sdf")
    dock_and_get_interactions_fn = partial(
        dock_and_get_interactions,
        protein_pdb_path=protein_pdb_path,
        autobox_ligand_path=autobox_ligand_path,
        residues_of_interest={
            "VAL629", "ARG715", "TRP718", "THR555", "CYS675", 
            # "LYS581", "GLN626", "LYS677", "SER698",
        },
        necessary_residues={"VAL629"},
        exhaustiveness=8,
        seed=1337,
        num_poses=5,
        temp_dir=Path(f"/home/minhtoo/molecule-generation/temp_dirs/moler_oracle_temp_{run_name}_{time_hash}/")
    )
    # TODO: add batched GNINA docking fn
    scoring_functions = [ScoringFunction(func=dock_and_get_interactions_fn, name="interactions", is_mol_func=True)]

    out_dir = Path(f"/home/minhtoo/molecule-generation/results/{run_name}_{time_hash}/")
    model_dir = Path("/home/minhtoo/molecule-generation/models/")
    with load_model_from_directory(model_dir) as model:
        infer_model = MoLeRModel(model, scaffold, debug=True)
        # TODO: possible to init with multiple smiles? need to refactor...
        opt = BasePSOptimizer.from_query(
            init_smiles=init_smiles,
            num_part=num_parts,
            num_swarms=1,
            inference_model=infer_model,
            scoring_functions=scoring_functions
        )
        opt.run(10, num_track=1_000_000, out_dir=out_dir)

    logger.success('best solutions')
    logger.success(opt.best_solutions)
    logger.success('best fitness history')
    logger.success(opt.best_fitness_history)


if __name__ == "__main__":
    optimise()