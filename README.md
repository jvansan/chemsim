# chemsim

A learning exercise in chemical similarity metrics.

Comparing RDKit similarity to self-implementations. All functions are design to handle RDKit Fingerprints.

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

Each version calculates the similairty of the same two fingerprints 1000 times and the mean times are:

```
**Tanimoto**
0.2785714285714286
0.2785714285714286
0.2785714285714286
0.2785714285714286
RDKit 2.559661865234375e-06
Numpy 0.006477154493331909
_popc 1.779031753540039e-05
Py 3.1856775283813474e-05

**Dice**
0.43575418994413406
0.43575418994413406
0.43575418994413406
RDKit 1.9543170928955076e-06
Numpy 0.006620363473892212
Py 3.6644935607910154e-05

**Cosine**
0.4357609900500941
0.4357609900500941
0.4357609900500941
RDKit 1.1813640594482422e-06
Numpy 0.011004742860794068
Py 3.5953760147094724e-05

**Soergel**
0.4357609900500941
0.4357609900500941
Numpy 0.006542242527008056
Py 3.6862850189208985e-05

**Manhattan**
0.4357609900500941
0.4357609900500941
Numpy 0.00639375376701355
Py 2.670145034790039e-05
```

#### Useful resources:

- http://www.dalkescientific.com/writings/diary/archive/2020/09/28/simple_fps_fingerprint_search.html
- https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics_OLCC_(2019)/6%3A_Molecular_Similarity/6.2%3A_Similarity_Coefficients
- https://www.sequentix.de/gelquest/help/distance_measures.htm
