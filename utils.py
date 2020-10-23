import numpy as np
import _popc
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from math import sqrt


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
    fp1 = DataStructs.BitVectToBinaryText(fp1)
    fp2 = DataStructs.BitVectToBinaryText(fp2)
    return _popc.lib.byte_tanimoto_256(fp1, fp2)


def popcount(x):
    return bin(x).count("1")


def py_tanimoto(fp1, fp2):
    """Using Python bitwise operations"""
    fp1 = int(fp1.ToBitString(), 2)
    fp2 = int(fp2.ToBitString(), 2)
    return popcount(fp1 & fp2) / popcount(fp1 | fp2)


def rdkit_dice(fp1, fp2):
    """Using builtin rdkit operations"""

    return DataStructs.DiceSimilarity(fp1, fp2)


def np_dice(fp1, fp2):
    """Using numpy logical operations"""
    fp1 = np.array(fp1)
    fp2 = np.array(fp2)

    return 2 * sum(np.logical_and(fp1, fp2)) / (sum(fp1) + sum(fp2))


def py_dice(fp1, fp2):
    """Using Python bitwise operations"""
    fp1 = int(fp1.ToBitString(), 2)
    fp2 = int(fp2.ToBitString(), 2)
    return 2 * popcount(fp1 & fp2) / (popcount(fp1) + popcount(fp2))


def py_soergel(fp1, fp2):
    """Using Python bitwise operations
    Soergel DISTANCE is the complement of Tanimoto Similarity.
    """
    fp1 = int(fp1.ToBitString(), 2)
    fp2 = int(fp2.ToBitString(), 2)
    a = popcount(fp1)
    b = popcount(fp2)
    c = popcount(fp1 & fp2)
    return (a + b - 2.0 * c) / (a + b - c)


def rdkit_cosine(fp1, fp2):
    """Using builtin rdkit operations"""

    return DataStructs.CosineSimilarity(fp1, fp2)


def np_cosine(fp1, fp2):
    """Using numpy logical operations"""
    fp1 = np.array(fp1)
    fp2 = np.array(fp2)

    return sum(np.logical_and(fp1, fp2)) / np.sqrt(sum(fp1) * sum(fp2))


def py_cosine(fp1, fp2):
    """Using Python bitwise operations"""

    fp1 = int(fp1.ToBitString(), 2)
    fp2 = int(fp2.ToBitString(), 2)

    return popcount(fp1 & fp2) / sqrt(1.0 * popcount(fp1) * popcount(fp2))
