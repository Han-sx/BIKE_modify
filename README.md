BIKE-1 - Additional implementation of "BIKE: Bit Flipping Key Encapsulation". 
-----------------------------------------------------------------------------
**【！before you read，the developing branch is the latest branch, please switch branches to view the codes！】**

This package is an "additional optimized" implementation of the Round-2 
variant of BIKE-1. 

BIKE is a submission to the Post-Quantum Cryptography 
Standardization project http://csrc.nist.gov/projects/post-quantum-cryptography.

The official BIKE website is: https://bikesuite.org

This package offers two constant time implementations of Round-2 BIKE-1.
- A portable implementation that requires libcrypto.a (e.g., of OpenSSL) for AES256 and SHA384.
- An x86_64 optimized implementation with AVX2 or AVX512.

The optimizations in this package are based on the papers: 
- Nir Drucker, Shay Gueron, "A Toolbox for Software Optimization of QC-MDPC Code-Based Cryptosystems." 
  Journal of Cryptographic Engineering, January 2019, 1–17 https://doi.org/10.1007/s13389-018-00200-4.

- Chou, T.: QcBits: Constant-Time Small-Key Code-Based Cryptography. In: 
  Gier-lichs, B., Poschmann, A.Y. (eds.) Cryptographic Hardware and 
  Embedded Systems– CHES 2016. pp. 280–300. Springer Berlin Heidelberg, 
  Berlin, Heidelberg (2016)

- Guimarães, Antonio, Diego F Aranha, and Edson Borin. 2019. 
  “Optimized Implementation of QC-MDPC Code-Based Cryptography.”
  Concurrency and Computation: Practice and Experience 31 (18): 
  e5089. https://doi.org/10.1002/cpe.5089.

The decoder (in decoder/decoder.c) algorithm is the Black-Gray decoder included
in the early submission of CAKE (due to N. Sandrier and R. Misoczki).

The analysis for the constant-time implementation is given in:
Drucker, Nir, Shay Gueron, and Dusan Kostic. 2019. 
“On Constant-Time QC-MDPC Decoding with Negligible Failure Rate.” 
Cryptology ePrint Archive, 2019. https://eprint.iacr.org/2019/1289.

The BG decoder was modified into BGF in
Drucker, Nir, Shay Gueron, and Dusan Kostic. 2019. “QC-MDPC decoders with several shades of gray.” 
Cryptology EPrint Archive, 2019.

The code is due to Nir Drucker, Shay Gueron, and Dusan Kostic.

BUILD
-----

To compile:
    make

Default is BIKE-1 (Round-2 variant) at Level-1 (64-bit quantum security) without AVX2/512 support (portable).
This compilation assumes that AES_NI, POPCNT, and PCLMULQDQ instructions are available, for other platforms use USE_OPENSSL=1.

Additional compilation flags:
 - USE_NIST_RAND - Using the RDBG of NIST and generate the KATs.
 - USE_OPENSSL   - Use OpenSSL for AES/SHA and GF2X multiplication. 
                   OpenSSL must be installed on the platform.
 - OPENSSL_DIR   - Set the path of the OpenSSL include/lib directories.
 - FIXED_SEED    - Using a fixed seed, for debug purposes.
 - RDTSC         - Measure time in cycles rather than in mseconds.
 - VERBOSE       - Add verbose (level:1-4 default:1).
 - NUM_OF_TESTS  - Set the number of tests to be run.
 - AVX2          - Compile with AVX2 support (to compile use GCC).
 - AVX512        - Compile with AVX512 support (to compile use GCC).
 - LEVEL         - Security level (1/3/5).
 - ASAN/TSAN - Enable the associated clang sanitizer
 
To clean:
    make clean.

Note: a "clean" is required prior to compilation with modified flags.

Example: 
    make AVX2=1 USE_OPENSSL=1 VERBOSE=2

The package was compiled and tested with gcc (version 4.8.0 or above) in 64-bit mode. 
Tests were run on a Linux (Ubuntu 16.04.3 LTS) OS. 
Compilation on other platforms may require some adjustments.
