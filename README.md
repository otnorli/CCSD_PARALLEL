RHF CCSD USING CLOSED SHELL GAUSSIAN TYPE ORBITALS FROM EMSL BASIS SET EXCHANGE - PROGRAM WRITTEN USING QTCREATOR

CCSDTQ program in "cc_general.cpp", it works with a different version of Hartree Fock and is not 100% completed jet, but uploaded still. Big program.

Big CCSD program. Very fast! Implemented brand new optimizations never before seen in the human world. 

Should be noted that it currently only works for even number of electrons, but we only need to add a few +1 and -1 around the code for it to work for odd number of electrons, no optimizations exclude odd number of electrons.

Eggstremely fast program, fastest in the world as far as I know and free to be downloaded by anyone! Feel free to use my code and make whatever changes you want. No need to ask permission. In fact please dont ask my permission.

ccsd_memory_optimized.cpp is the most optimal code. Rest of CCSD classes are stored versions of earlier and less effective code, in both speed and memory terms.

Program not 100% complete, still in the works. 

Program is probably hard to get to work right now, as input can bug in meny meny ways. But it will be much easier soon. User friendliness will be good.

Program uses Armadillo and MPI. No OpenMP or any of that. Please supply blas and lapack also somehow. Input is given in a file you name "INCAR". 

Here supply convergence criteria, Bohr og angstrom, atom positions and atom names etc. Description and examples comming soon.
