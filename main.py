import time
import utils

# cids = [    54454,  # Simvastatin (Zocor)
#             54687,  # Pravastatin (Pravachol)
#             60823,  # Atorvastatin (Lipitor)
#            446155,  # Fluvastatin (Lescol)
#            446157,  # Rosuvastatin (Crestor)
#           5282452,  # Pitavastatin (Livalo)
#          97938126 ] # Lovastatin (Altoprev)

COMPOUNDS = [
    "CCC(C)(C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
    "CC[C@H](C)C(=O)O[C@H]1C[C@@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@H](C[C@H](CC(=O)O)O)O)O",
    "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
    "CC(C)N1C2=CC=CC=C2C(=C1/C=C/[C@H](C[C@H](CC(=O)O)O)O)C3=CC=C(C=C3)F",
    "CC(C)C1=NC(=NC(=C1/C=C/[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
    "C1CC1C2=NC3=CC=CC=C3C(=C2/C=C/[C@H](C[C@H](CC(=O)O)O)O)C4=CC=C(C=C4)F",
    "CC[C@H](C)C(=O)O[C@@H]1C[C@@H](C[C@@H]2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
]


def mean_timing(fn, max_iter=1000):
    times = []
    for _ in range(max_iter):
        start = time.time()
        fn()
        end = time.time()
        times.append(end - start)
    return sum(times) / len(times)


def main():
    fps = [utils.smi2mfp4(c) for c in COMPOUNDS]

    rdk_time = mean_timing(lambda: utils.rdkit_tanimoto(fps[0], fps[1]))
    np_time = mean_timing(lambda: utils.np_tanimoto(fps[0], fps[1]))
    popc_time = mean_timing(lambda: utils.popc_tanimoto(fps[0], fps[1]))
    py_time = mean_timing(lambda: utils.py_tanimoto(fps[0], fps[1]))
    print("RDKit", rdk_time)
    print("Numpy", np_time)
    print("_popc", popc_time)
    print("Py", py_time)


if __name__ == "__main__":
    main()
