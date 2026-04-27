import argparse
import io
import json
import logging
import os
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Optional, Type, TypeVar, Union

import MDAnalysis as mda
import yaml
from openbabel import pybel
from pydantic import BaseModel as _BaseModel
from pydantic import validator
from rdkit import Chem, RDLogger

T = TypeVar("T")
PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def configure_logging(debug=0) -> None:
    """Set up logging."""
    logger_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        level=logger_level,
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )


def _resolve_path_exists(value: Optional[Path]) -> Optional[Path]:
    """Check if a path exists (implements path_validator)."""
    if value is None:
        return None
    p = value.resolve()
    if not p.exists():
        raise FileNotFoundError(p)
    return p


def _resolve_mkdir(value: Path) -> Path:
    """Create a directory if it does not exist (implements mkdir_validator)."""
    p = value.resolve()
    p.mkdir(exist_ok=False, parents=True)
    return p.absolute()


def path_validator(field: str) -> classmethod:
    """Pydantic validator to check if a path exists."""
    decorator = validator(field, allow_reuse=True)
    _validator = decorator(_resolve_path_exists)
    return _validator


def mkdir_validator(field: str) -> classmethod:
    """Pydantic validator to create a directory if it does not exist."""
    decorator = validator(field, allow_reuse=True)
    _validator = decorator(_resolve_mkdir)
    return _validator


class BaseModel(_BaseModel):
    """An interface to add JSON/YAML serialization to Pydantic models"""

    def write_json(self, path: PathLike) -> None:
        """Write the model to a JSON file.

        Parameters
        ----------
        path : str
            The path to the JSON file.
        """
        with open(path, "w") as fp:
            json.dump(self.dict(), fp, indent=2)

    @classmethod
    def from_json(cls: Type[T], path: PathLike) -> T:
        """Load the model from a JSON file.

        Parameters
        ----------
        path : str
            The path to the JSON file.

        Returns
        -------
        T
            A specific BaseModel instance.
        """
        with open(path, "r") as fp:
            data = json.load(fp)
        return cls(**data)

    def write_yaml(self, path: PathLike) -> None:
        """Write the model to a YAML file.

        Parameters
        ----------
        path : str
            The path to the YAML file.
        """
        with open(path, mode="w") as fp:
            yaml.dump(json.loads(self.json()), fp, indent=4, sort_keys=False)

    @classmethod
    def from_yaml(cls: Type[T], path: PathLike) -> T:
        """Load the model from a YAML file.

        Parameters
        ----------
        path : PathLike
            The path to the YAML file.

        Returns
        -------
        T
            A specific BaseModel instance.
        """
        with open(path) as fp:
            raw_data = yaml.safe_load(fp)
        return cls(**raw_data)  # type: ignore

    @classmethod
    def from_args(cls: Type[T], args: argparse.Namespace) -> T:
        """Load the model from an argparse result

        Parameters
        ----------
        args: argparse.Namespace
            The argparse argument set.

        Returns
        -------
        T
            A specific BaseModel instance.
        """
        v = {}
        adict = vars(args)
        raw = {"dock": v}
        for key in [
            "name",
            "receptor_pdb",
            "drug_dbs",
            "diffdock_dir",
            "output_dir",
            "top_n",
            "batch_size",
        ]:
            if key in args:
                v[key] = adict[key]
        return cls(**raw)  # type: ignore


three_to_one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",  # MSE this is almost the same AA as MET. The sulfur is just replaced by Selen
    "PHE": "F",
    "PRO": "P",
    "PYL": "O",
    "SER": "S",
    "SEC": "U",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "XAA": "X",
    "XLE": "J",
}


# def build_logger(debug=0):
#     logger_level = logging.DEBUG if debug else logging.INFO
#     logging.basicConfig(level=logger_level, format="%(asctime)s %(message)s")
#     logger = logging.getLogger(__name__)
#     return logger


def run_and_save(cmd, cwd=None, output_file=None):
    print(cmd, file=output_file)
    cmd = cmd.split()
    process = subprocess.Popen(
        cmd,
        shell=False,
        cwd=cwd,
        stdout=output_file,
        stderr=subprocess.STDOUT,
    )
    rc = process.wait()
    if rc != 0:
        print(f"Failure running {cmd}", file=sys.stderr)
        sys.exit(1)
    return process

def run_list_and_save(cmd, cwd=None, output_file=None):
    print(cmd, file=output_file)
    process = subprocess.Popen(
        cmd,
        shell=False,
        cwd=cwd,
        stdout=output_file,
        stderr=subprocess.STDOUT,
    )
    rc = process.wait()
    if rc != 0:
        print(f"Failure running {cmd}", file=sys.stderr)
        sys.exit(1)
    return process


def pdb2seq(pdb_file):
    mda_u = mda.Universe(pdb_file)
    protein = mda_u.select_atoms("protein")
    seq = ""
    for res in protein.residues:
        if res.resname in three_to_one:
            seq += three_to_one[res.resname]
        else:
            seq += "-"
    return seq


def get_pdblabel(pdb_file) -> str:
    return os.path.basename(pdb_file)[:-4]


def comb_pdb(prot_pdb, lig_pdb, comp_pdb=None):
    if comp_pdb is None:
        comp_pdb = f"{os.path.dirname(lig_pdb)}/{get_pdblabel(prot_pdb)}_{get_pdblabel(lig_pdb)}.pdb"
    if os.path.exists(comp_pdb):
        os.remove(comp_pdb)

    prot_u = mda.Universe(prot_pdb)
    lig_u = mda.Universe(lig_pdb)
    merged = mda.Merge(prot_u.atoms, lig_u.atoms)
    try:
        merged.atoms.write(comp_pdb)
        return comp_pdb
    except Exception as e:
        logging.info(f"Skipping {lig_pdb} for position error {e}. ")
        return None


def sdf2pdb(sdf_file, pdb_file=None):
    if pdb_file is None:
        pdb_file = sdf_file[:-3] + "pdb"

    if os.path.exists(pdb_file):
        os.remove(pdb_file)

    try:
        mol = next(pybel.readfile("sdf", sdf_file))
        mol.write("pdb", pdb_file)
    except Exception as err:
        print(
            f"Failure (output is empty) converting f{sdf_file} to f{pdb_file} err={err}"
        )
        sys.exit(1)

    return pdb_file


def clean_pdb(pdb_file: str, output_pdb: str) -> str:
    if pdb_file.endswith("cif.gz") or pdb_file.endswith("cif"):
        logger.debug(
            f"Found cif format for input target protein. Converting it to pdb. "
        )
        mol = next(pybel.readfile("cif", pdb_file))

        tempdir = tempfile.TemporaryDirectory()
        temp_pdb = (
            f"{tempdir.name}/{os.path.basename(pdb_file.replace('.gz', ''))[:-4]}.pdb"
        )
        mol.write("pdb", temp_pdb)
    else:
        temp_pdb = pdb_file

    mda_u = mda.Universe(temp_pdb)
    if mda_u.select_atoms("nucleic").n_atoms > 0:
        logger.warning(
            f"{os.path.basename(pdb_file)} contains nucleic acids (DNA/RNA), which will be removed due to incompatibility with the docking program. Proceed with caution. "
        )

    protein = mda_u.select_atoms("protein")
    if protein.n_residues > 1024:
        logger.warning(
            f"{os.path.basename(pdb_file)} has large number of residues (>1024). Proceed with caution. "
        )
    protein.write(output_pdb)
    return output_pdb


#
# verify that a smiles string is a smiles string.
#


def validate_smiles(smiles_str):
    is_valid = False
    RDLogger.DisableLog("rdApp.error")
    if Chem.MolFromSmiles(smiles_str) is not None:
        is_valid = True
    RDLogger.EnableLog("rdApp.error")

    return is_valid


def cal_cnn_aff(pdb_file, sdf_file, gnina_exe="gnina", log_handle=None):
    tempdir = tempfile.TemporaryDirectory()

    output_sdf = f"{tempdir.name}/{os.path.basename(sdf_file)}"
    # Use --no_gpu for CPU-only mode to avoid CUDA compatibility issues
    cmd = (
        f"{gnina_exe} --no_gpu --minimize --scoring vinardo "
        f"-r {pdb_file} -l {sdf_file} "
        f"--autobox_ligand {sdf_file} "
        f"--autobox_add 2 -o {output_sdf}"
    )
    run_and_save(cmd, output_file=log_handle)

    try:
        mol = next(pybel.readfile("sdf", output_sdf))
        return mol
    except Exception as e:
        print(f"Failed gnina run on {sdf_file} with {e}. ")
        return None


# def prep_receptor(receptor_file: str) -> str:
#     if receptor_file.endswith("cif.gz") or receptor_file.endswith("cif"):
#         mol = next(pybel.readfile("cif", receptor_file))

#         tempdir = tempfile.TemporaryDirectory()
#         temp_pdb = f"{tempdir.name}/{os.path.basename(receptor_file.replace('.gz', ''))[:-3]}.pdb"
#         mol.write("pdb", temp_pdb)
#         return temp_pdb
#     else:
#         return os.path.abspath(receptor_file)
