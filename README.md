# chemsim

A learning exercise in chemical similarity metrics.

Comparing RDKit similarity to self-implementations

#### Current results

Each version calculates the tanimoto of the same two fingerprints 1000 times and the mean times are:

```
RDKit 2.2430419921875e-06
Numpy 0.013622828722000123
_popc 1.1498451232910157e-05
Py 0.00011100435256958008
```
