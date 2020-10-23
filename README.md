# chemsim

A learning exercise in chemical similarity metrics.

Comparing RDKit similarity to self-implementations

#### Current results

Each version calculates the similairty of the same two fingerprints 1000 times and the mean times are:

```
**Tanimoto**
0.2785714285714286
0.2785714285714286
0.2785714285714286
0.2785714285714286
RDKit 1.935243606567383e-06
Numpy 0.01377569031715393
_popc 1.678895950317383e-05
Py 3.043198585510254e-05

**Dice**
0.43575418994413406
0.43575418994413406
0.43575418994413406
RDKit 2.1321773529052733e-06
Numpy 0.010457364082336426
Py 3.805208206176758e-05

**Cosine**
0.4357609900500941
0.4357609900500941
0.4357609900500941
RDKit 1.2481212615966798e-06
Numpy 0.010791512012481689
Py 3.4531831741333005e-05
```
