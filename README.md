Quantum Chemistry Package from Computational Physics Group at University of Oslo.

Available Methods:

- UHF, Unstricted Hartree-Fock
- RHF, Restricted Hartree Fock
- CCSD, Coupled Cluster Singles and Doubles, with spin restriction
- CCSDT, Coupled Cluster Singles, Doubles and Triples, with spin restriction

and CCSDT-n Models:
- CCSDT-1a
- CCSDT-1b
- CCSDT-2
- CCSDT-3
- CCSDT-4

Implemented with Cartesian Gaussian Type Orbitals, Cartesian GTOs.



Please see file "INCAR". This is input file for program.



Description of classes in alphabetical order:

- cc_general.cpp:
Unfinished CCSDTQ Method.

- ccsd_even.cpp:
Old version of CCSD program using Compact Storage of arrays.

- ccsd_memory_optimized.cpp
Parallel working memory distributed highly optimized superhighspeed CCSD program. 

- ccsd_non_iterative_part.cpp
Class used to optimize ccsD_memory_optimized further. Maps out MPI.

- ccsdt.cpp
Calculates CCSDT energies.

- coupled_cluster.cpp:
First implementation of CCSD.

- fill_alpha.cpp:
Contains numbers from EMSL for a select number of basis sets.

- hartree_fock_solver.cpp:
Solves the Hartree Fock Equations for RHF and UHF.

- hartree_integrals:
Solves the Overlap, Kinetic, Nuclei-Electron and Electron-Electron integrals using Cartesian GTOs.

- initializer.cpp:
Makes use of input to run correct methods with atoms in correct geometries.

- main.cpp:
Reads input file and sends information to initializer.

- matrix_size_setter:
Sets up sizes of arrays we will need.



Additional Comments:

- CCSD Program is really a very high performance implementation. Still somewhat sub optimal (relative to perfect) parallel performance. Also some scaling problems with armadillo in combination with MPI. Utilizes newly discovered compact storage to effectively dodge all uneccasary calculations and utilize symmetries were they are effective. 

- For further proposed optimizations please see MA thesis "Coupled Cluster Studies in Computational Chemistry" from the Computational Physics Group at the University of Oslo. Anyone can copy any parts of this code.

- Makes use of armadillo and MPI. No other external libraries required.



User guide:

1) Install armadillo and MPI.

2) Ensure armadillo and MPI is linked properly. Edit .h files if needed. 

3) Ensure armadillo is installed properly to make use of external math library, like Intel MKL or BLAS, if you want high performance.

4) Input given in input file, which must be named INCAR. Please see uploaded example. 

5) If the combination of atoms you want are not currently available, try a different basis set or update fill_alpha and matrix_size_setter classes with new numbers from EMSL Basis Set Exchange (Freely available on the internet, google it).


