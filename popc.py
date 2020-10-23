# Code from http://www.dalkescientific.com/writings/diary/archive/2020/09/28/simple_fps_fingerprint_search.html
# I call this file "popc.py". It should work for gcc and clang.
from cffi import FFI

ffibuilder = FFI()

# Create a Python module which can be imported via "import _popc".
# It will contain the single function byte_tanimoto_256() which
# expects two byte strings of length 256 bytes - exactly!
ffibuilder.set_source(
    "_popc",
    r"""

#include <stdint.h>

static double
byte_tanimoto_256(const unsigned char *fp1, const unsigned char *fp2) {
    int num_words = 2048 / 64;
    int union_popcount = 0, intersect_popcount = 0;

    /* Interpret as 64-bit integers and assume possible mis-alignment is okay. */
    uint64_t *fp1_64 = (uint64_t *) fp1, *fp2_64 = (uint64_t *) fp2;

    for (int i=0; i<num_words; i++) {
        intersect_popcount += __builtin_popcountll(fp1_64[i] & fp2_64[i]);
        union_popcount += __builtin_popcountll(fp1_64[i] | fp2_64[i]);
    }
    if (union_popcount == 0) {
        return 0.0;
    }
    return ((double) intersect_popcount) / union_popcount;
}

""",
    # Tell the compiler to always expect the POPCNT instruction will be available.
    extra_compile_args=["-mpopcnt"],
)

# Tell cffi to export the above function for use by Python.
ffibuilder.cdef(
    """
double byte_tanimoto_256(unsigned char *fp1, unsigned char *fp2);
"""
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
