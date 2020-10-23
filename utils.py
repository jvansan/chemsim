import numpy as np
import _popc
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def smi2mfp4(smi):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, 4, nBits=2048)


def rdkit_tanimoto(fp1, fp2):
    """Using builtin rdkit operations"""

    return DataStructs.TanimotoSimilarity(fp1, fp2)


def np_tanimoto(fp1, fp2):
    """Using numpy logical operations"""
    fp1 = np.array(fp1)
    fp2 = np.array(fp2)

    return sum(np.logical_and(fp1, fp2)) / sum(np.logical_or(fp1, fp2))


def popc_tanimoto(fp1, fp2):
    """Using C wrapped bitwise operations
    Code from: http://www.dalkescientific.com/writings/diary/archive/2020/09/28/simple_fps_fingerprint_search.html
    Requires you to run `python popc.py` to compile
    """
    return _popc.lib.byte_tanimoto_256(fp1.ToBinary(), fp2.ToBinary())


def popcount(x):
    return bin(x).count("1")


def py_tanimoto(fp1, fp2):
    """Using Python bitwise operations"""
    fp1 = int(fp1.ToBitString())
    fp2 = int(fp2.ToBitString())
    return popcount(fp1 & fp2) / popcount(fp1 | fp2)

