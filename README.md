# polar-decomposition-3x3
Implementation of polar decomposition for 3x3 matrices

This project is a C++ implementation of polar decomposition for 3x3 matrices.

It is based on the paper:

N. J. Higham and V. Noferini. **"An algorithm to compute the polar decomposition of a 3x3 matrix"**. Numer. Algorithms, 73(2):349-369, 2016.

The paper is available at:

- http://eprints.ma.man.ac.uk/2352/01/covered/MIMS_ep2015_66.pdf

It is also now available at:

- http://link.springer.com/article/10.1007%2Fs11075-016-0098-7

This C++ implementation is also extensively based on the Matlab implementation of this algorithm available at:

- https://github.com/higham/polar-decomp-3by3

The implementation is a header-only library, so to use it in any given project, simply:

- copy the header files from the `include/` folder into the project
- include the `"polar_decomposition_3x3.h"` file in a source file
- call the `polar_decomposition()` function.

The implementation is templated so computation can be performed either in single or double precision.  It expects matrices in column-major representation.  A sample test file is provided to demonstrate how to include and call the library.

It is worth noting that Matlab uses a `(row,column)` matrix indexing, but this implementation uses `(column,row)` indexing.  Also, Matlab is 1-index based, but this implementation is 0-index based.

This implementation has very specific goals that drive implementation choices:

- It tries to avoid dependencies to third-party libraries.  Therefore, it reimplements basic operations (matrix operations such as multiplication, transposition, etc.).  These are straight-forward to implement and are kept separate to the actual algorithm so that using an actual linear algebra library would make the implementation more straight-forward.
- It tries to be as efficient as possible.  Therefore it potentially combines multiple operations into one (for instance, combine matrix transposition with multiplication) to minimize runtime cost.  While this might reduce code simplicity, we favor runtime efficiency while trying to making those optimizations as easy to read as possible.

It is also worth noting that this implementation relies on two major sources:

- The algorithm as described in the paper
- The algorithm as implemented in the Matlab implementation.

This implementation tries to highlight its references to both the paper and the source code.
