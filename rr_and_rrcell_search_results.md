# Files involving RR and RRcell in Macaulay2/d and Macaulay2/e directories

## Summary

This document lists all files found in the `M2/Macaulay2/d` and `M2/Macaulay2/e` directories that contain references to `RR` (real numbers) and `RRcell` (real number cell type) patterns.

## Files in M2/Macaulay2/d directory

### Files with RR patterns:

1. **basic.d** - Contains hash function for RRcell
2. **util.d** - Core RR/RRcell conversion functions
3. **texmacs.d** - Error handling with stderr
4. **tokens.d** - Error message printing functions
5. **stdio0.d** - STDERR constant definition
6. **stdio.d** - String conversion and output functions for RR
7. **scclib.c** - Error output to stderr
8. **regex.dd** - Exception handling with stderr output
9. **python.d** - Python interface functions for RR/RRcell
10. **pthread.d** - Thread management with error handling
11. **parse.d** - Parser constants and real number code definitions
12. **parser.d** - RR parsing functions
13. **mysqldummy.d** - MySQL error handling functions
14. **memdebug.c** - Memory debugging with stderr output
15. **lex.d** - Lexer with RR type codes
16. **interface2.d** - Interface functions with RR argument validation
17. **interface.dd** - Interface functions with RR conversions
18. **actors2.dd** - Mathematical operations on RRcell
19. **actors3.d** - Extensive mathematical functions for RRcell (trigonometric, logarithmic, special functions)
20. **actors4.d** - Formatting and conversion functions for RRcell
21. **actors5.d** - Additional mathematical functions for RRcell
22. **evaluate.d** - Expression evaluation with RRcell

### Files with RRcell patterns:

1. **actors2.dd** - Basic RRcell operations and type checking
2. **util.d** - RRcell creation and conversion functions
3. **actors3.d** - Mathematical operations on RRcell (sin, cos, exp, log, Gamma, Bessel functions, etc.)
4. **python.d** - Python interface for RRcell
5. **parse.d** - Parser type definitions including RRcell
6. **interface2.d** - Interface functions with RRcell handling
7. **basic.d** - Hash functions for RRcell
8. **actors5.d** - Additional RRcell operations (factorial, external string conversion)
9. **evaluate.d** - RRcell expression evaluation
10. **actors4.d** - RRcell formatting, precision, and conversion functions
11. **interface.dd** - RRcell to engine interface conversions

## Files in M2/Macaulay2/e directory

### Files with RR patterns:

1. **Makefile.files** - Build system references to aring-RR and aring-RRR
2. **TODO-rings-matrices** - Documentation about RR, RRR ring implementations
3. **aring-gf-flint.hpp** - BigReal conversion functions
4. **aring-qq-flint.hpp** - BigReal conversion functions
5. **debug.hpp** - Debug function for RRR
6. **SLP-defs.hpp** - Homotopy algorithm definitions for RR/RRR
7. **dmat-ffpack.cpp** - Linear algebra error output
8. **dmat-lu-inplace.hpp** - LU decomposition for RR/RRR matrices
9. **finalize.cpp** - Object finalization with stderr output
10. **matrix.hpp** - Matrix operations for RRR/CCC
11. **eigen.cpp** - Eigenvalue/SVD computations for RR matrices
12. **interface/aring.h** - Ring interface functions
13. **skewpoly.cpp** - Skew polynomial ring initialization
14. **interface/groebner.h** - Gröbner basis computations with RR support
15. **LLL.hpp** - LLL algorithm documentation
16. **Makefile.in** - Build configuration
17. **CMakeLists.txt** - CMake build system with RR-related targets
18. **interface/aring.cpp** - Ring arithmetic implementations

### Files with RRcell patterns:

No files in the `M2/Macaulay2/e` directory contain `RRcell` patterns. This suggests that `RRcell` is primarily a data structure used in the D language frontend (`d` directory), while the engine (`e` directory) works with different representations of real numbers.

## Key Observations

1. **RRcell is frontend-specific**: The `RRcell` type appears to be used exclusively in the D language frontend code in the `d` directory for representing real number values in expressions.

2. **RR in engine**: The `e` directory (engine) uses `RR`, `RRR` and related types for actual computational arithmetic and linear algebra operations.

3. **Mathematical functions**: Most mathematical functions (trigonometric, logarithmic, special functions) are implemented with `RRcell` handling in the `actors3.d` file.

4. **Interface layer**: The `interface.dd` and `interface2.d` files handle conversion between the frontend `RRcell` representation and the engine's internal representations.

5. **Build system integration**: The engine build system (`CMakeLists.txt`, `Makefile.files`) includes specific targets for RR arithmetic implementations.