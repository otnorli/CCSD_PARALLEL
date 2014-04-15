RHF CCSD USING CLOSED SHELL GAUSSIAN TYPE ORBITALS FROM EMSL BASIS SET EXCHANGE

Big CCSD program. Very fast! Implemented brand new optimizations never before seen in the human world. 

Should be noted that it currently only works for even number of electrons, but we only need to add a few +1 and -1 around the code for it to work for odd number of electrons, no optimizations exclude odd number of electrons.

Eggstremely fast program, fastest in the world as far as I know and free to be downloaded by anyone! Feel free to use my code and make whatever changes you want. No need to ask permission. In fact please dont ask my permission.

ccsd_memory_optimized.cpp is the most optimal code. Rest of CCSD classes are stored versions of earlier and less effective code, in both speed and memory terms.

Program not 100% complete, still in the works. Simple minded version of CCSDTQ soon to be uploaded, also to be freely available.

Program uses Armadillo and MPI. No OpenMP or any of that. Please supply blas and lapack also. Input is given in a file you name "INCAR". In this file you write:

\#ATOMS BEGIN

then list a bunch of atoms, like 

O 0 0 0

H 1 0 0
 etc

and then \#ATOMS END

also add a few lines about stuff like
Basis_Set STO-3G
Method CCSD
use_angstrom true
convergance_criteria -8.0

If you want to change basis set you can use for example 4-31G, STO-3G, 6-311-2d2p or 6-311ss. Not all atoms implemented for all these basis sets jet, but soon.

The only important thing is that you have \#ATOMS START and \#ATOMS END and in the middle you have the chemical symbol for the grunnstoff you wish to study, O for oxygen, H for hydrogen etc, and then you have the coordinates x y z with one space in between them. Or else you will get bug. Can give input in main.cpp also if you edit some code. 

Output will be printed if print_stuffies is set to true.
Ångstrøm is used if use_angstrom is set to true
relax_pos does not work jet, set it false
convergance_criteria is given as the power of 10, for example -8.0 is 10^(-8.0) convergance criteria
Method can be CCSD or HF currently
Basis_Set are listed above
