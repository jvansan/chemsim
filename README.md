# chemsim

A learning exercise in chemical similarity metrics.

Comparing RDKit similarity to self-implementations. All functions are designed to handle RDKit Fingerprints.

Metrics implemented:

1. Tanimoto
2. Dice
3. Cosine
4. Soergel\*
5. Euclidean\*
6. Manhattan\*

\* Distances converted to similarity as per [Equation 1 in ref.](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3)

**Equation 1**

```
similarity = 1 / (1 + distance)
```

#### Current results

Each version calculates the similairty of the same set of fingerprints 1000 times and the mean times are:

```
**Tanimoto**
0.2785714285714286
0.2785714285714286
0.2785714285714286
0.2785714285714286
RDKit 2.0730836050851004e-06
Numpy 0.006472579956054688
_popc 1.7363718577793666e-05
Py 3.146038736615862e-05

**Dice**
0.43575418994413406
0.43575418994413406
0.43575418994413406
RDKit 1.8496853964669365e-06
Numpy 0.006422779866627285
Py 3.7190130778721404e-05

**Cosine**
0.4357609900500941
0.4357609900500941
0.4357609900500941
RDKit 1.0515621730259486e-06
Numpy 0.011113291604178292
Py 3.9301974432809014e-05

**Soergel**
0.5809128630705394
0.5809128630705394
Numpy 0.006626066752842495
Py 3.655951363699777e-05

**Manhattan**
0.00980392156862745
0.00980392156862745
Numpy 0.006840519768851144
Py 2.8515747615269253e-05

**Euclidean**
0.09049875621120891
0.09049875621120891
Numpy 0.006743010589054653
Py 2.8675828661237443e-05
```

#### Useful resources:

- http://www.dalkescientific.com/writings/diary/archive/2020/09/28/simple_fps_fingerprint_search.html
- https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics_OLCC_(2019)/6%3A_Molecular_Similarity/6.2%3A_Similarity_Coefficients
- https://www.sequentix.de/gelquest/help/distance_measures.htm
