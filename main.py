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


def main():
    fps = [utils.smi2mfp4(c) for c in COMPOUNDS]

    def mean_timing(fn, max_iter=1000):
        times = []
        for _ in range(max_iter):
            for fp in fps:
                start = time.time()
                fn(fps[0], fp)
                end = time.time()
                times.append(end - start)
        return sum(times) / len(times)

    # print("**Tanimoto**")
    # print(utils.rdkit_tanimoto(fps[0], fps[1]))
    # print(utils.np_tanimoto(fps[0], fps[1]))
    # print(utils.popc_tanimoto(fps[0], fps[1]))
    # print(utils.py_tanimoto(fps[0], fps[1]))
    # rdk_time = mean_timing(utils.rdkit_tanimoto)
    # np_time = mean_timing(utils.np_tanimoto)
    # popc_time = mean_timing(utils.popc_tanimoto)
    # py_time = mean_timing(utils.py_tanimoto)
    # print("RDKit", rdk_time)
    # print("Numpy", np_time)
    # print("_popc", popc_time)
    # print("Py", py_time)

    # print("**Dice**")
    # print(utils.rdkit_dice(fps[0], fps[1]))
    # print(utils.np_dice(fps[0], fps[1]))
    # print(utils.py_dice(fps[0], fps[1]))
    # rdk_time = mean_timing(utils.rdkit_dice)
    # np_time = mean_timing(utils.np_dice)
    # py_time = mean_timing(utils.py_dice)
    # print("RDKit", rdk_time)
    # print("Numpy", np_time)
    # print("Py", py_time)

    # print("**Cosine**")
    # print(utils.rdkit_cosine(fps[0], fps[1]))
    # print(utils.np_cosine(fps[0], fps[1]))
    # print(utils.py_cosine(fps[0], fps[1]))
    # rdk_time = mean_timing(utils.rdkit_cosine)
    # np_time = mean_timing(utils.np_cosine)
    # py_time = mean_timing(utils.py_cosine)
    # print("RDKit", rdk_time)
    # print("Numpy", np_time)
    # print("Py", py_time)

    print("**Soergel**")
    print(utils.np_soergel_sim(fps[0], fps[1]))
    print(utils.py_soergel_sim(fps[0], fps[1]))
    np_time = mean_timing(utils.np_soergel_sim)
    py_time = mean_timing(utils.py_soergel_sim)
    print("Numpy", np_time)
    print("Py", py_time)

    # print("**Manhattan**")
    # print(utils.np_manhattan_sim(fps[0], fps[1]))
    # print(utils.py_manhattan_sim(fps[0], fps[1]))
    # np_time = mean_timing(utils.np_manhattan_sim)
    # py_time = mean_timing(utils.py_manhattan_sim)
    # print("Numpy", np_time)
    # print("Py", py_time)

    # print("**Euclidean**")
    # print(utils.np_euclidean_sim(fps[0], fps[1]))
    # print(utils.py_euclidean_sim(fps[0], fps[1]))
    # np_time = mean_timing(utils.np_euclidean_sim)
    # py_time = mean_timing(utils.py_euclidean_sim)
    # print("Numpy", np_time)
    # print("Py", py_time)


if __name__ == "__main__":
    main()
