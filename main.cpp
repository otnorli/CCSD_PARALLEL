#include <iostream>
#define ARMA_NO_DEBUG
#include "../../home/ole/Desktop/include/armadillo"
#include <initializer.h>
#include <time.h>
#include <string>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace arma;
using namespace std;

/*
  This is a program applying the Restricted Closed Shell Hartree-Fock method and CCSD to determine an approximation to the
  wavefunction and the ground state energy.
  This program also has DIIS implemented for speedup in RHF and this is always used (unless you change the code somewhat).
  This can sometimes give bug, soon try: and catch: will fix this but not implemented jet

  Assumptions in HF:
  - Born-Opphenheimer (?)
  - Mean field
  - Relativistic effects neglected
  - The variational solution is assumed to be a linaer combination of a finite number of basis functions, which is chosen to be orthogonal
  - The combination of Gaussian Type Orbitals (GTOs) which is used is supposed to in combination resembel Slater orbitals

  Syntax:
  - C++
  - Library package Armadillo
  - Some oldschool C

  Comment language:
  - Mostly english
  - Some norwegian

  Inner workings:
  - Filling up alpha automaticly, or manually if you want, add it manually in fill_alpha.cpp. Just do the same as is done
  - Implemented meny basis sets, should be reasonably easy to understand how it works from the if tests
  - See hartree_fock_solver.cpp for the algorithm and hartree_integrals.cpp for the mathematics
  - ccsd_memory_optimized is the complete CCSD method, all in one class for easy and fast implementation

  Code is big and probably hard to understand. Soon there will be a paper or text describing the code, it is very easy its just a lot of
  symmetries and some communication optimization that makes it appear complex

  Brukermanual:
  - Ikke kjør så veldig rare systemer, foreløpig er ikke alt av basis set implementert for alle atomer
  - Bruk gjerne STO-3G eller 4-31G, de funker for det meste

  Input to give:
  - Method
  - Basis set
  - Nuclei charges and positions in a.u.
  - Extras

  Input given in external file INCAR, or change code modestly here in main file (uncomment and comment some lines) and you can give input here in main.cpp

  No need to supply things like how meny atoms, bond distances, or any of that. Just uncomment a method, uncomment a basis set,
  write some charges in Z and their positions in R. Very easy, some code has been implemented to simplify the input, this might
  seem confusing at first sight, but makes very easy and simple input.

  Current small bugs:
  - DIIS buggy sometimes, "Failed to solve" error = turn of DIIS in HF :O

  Program written by Ole Norli
*/

struct one_atom
{
    int CHARGE;
    double POS1;
    double POS2;
    double POS3;
    int CHECK_NEXT;
};

int return_charge(string Atom_Type);

char *find_and_write(char *fileBuffer, char *substring)
{
        char *ptr;
        int len = 0;

        if ( (ptr = strstr(fileBuffer, substring)) == NULL)
                return NULL; /* Error checking hvis return null etc */

        ptr += strlen(substring) + 1;

        // Sjekk lengde
        char *_ptr;
        for (_ptr = ptr; *_ptr != '\n'; _ptr++, len++)
                ;

        char *retval = (char *)malloc(len);
        memcpy(retval, ptr, len);
        *(retval + len) = '\0';
        return retval;
}

char *find_start_of_atoms(char *fileBuffer)
{
    // NEED "#ATOMS START" as a line in the INCAR file
    char *substring = "#ATOMS START";
    char *ptr;
    int len = 0;

    if ( (ptr = strstr(fileBuffer, substring)) == NULL)
            return NULL; /* Error checking hvis return null etc */

    ptr += strlen(substring) + 1;

    // Sjekk lengde
    char *_ptr;
    for (_ptr = ptr; *_ptr != '\n'; _ptr++, len++)
            ;
    return ptr;
}

void Fill_One_Atom(char**fileBuffer, one_atom *atom1)
{
    // BRUTE FORCE WAY of getting all atoms out, we take one line from input, seperate the values out
    // then delete the values we have seperated and return what remains in fileBuffer.
    // We do it this way because if we for example look for a certain atom , like "H", there may be more than one
    // of these atoms, that means we get problems having the right positions for the right atoms.

    // We also need to NOT GIVE ERROR if someone presses space additional times afterwards etc
    // However input does require a # sign at the end, like #ATOMS END or #END ATOMS or something
    // Also require this at start, an #ATOMS START

    char *ptr = *fileBuffer;
    int len = 0;
    double pos1, pos2, pos3;
    char *_ptr;

    // Sjekk lengde av første input ting, bokstavene som definerer hvem atom vi har
    for (_ptr = ptr; *_ptr != 0x20; _ptr++, len++)
            ;

    char *value = (char *)malloc(len);
    memcpy(value, *fileBuffer, len);
    *(value + len) = '\0';

    int atom_charge = return_charge(value);

    // We remove one and one part of the input file and add this position
    _ptr++;
    ptr = _ptr;
    len = 0;
    for (_ptr = ptr; *_ptr != 0x20; _ptr++, len++) // (0x20) = a (space), or (' ')
            ;

    value = (char *)realloc(value, len);
    memcpy(value, ptr, len);
    *(value + len) = '\0';
    pos1 = (double) atof(value);

    // Next
    _ptr++;
    ptr = _ptr;
    len = 0;
    for (_ptr = ptr; *_ptr != 0x20; _ptr++, len++)
            ;

    value = (char *)realloc(value, len);
    memcpy(value, ptr, len);
    *(value + len) = '\0';
    pos2 = (double) atof(value);

    // Next
    _ptr++;
    ptr = _ptr;
    len = 0;
    for (_ptr = ptr; *_ptr != 0x20; _ptr++, len++)
            ;

    value = (char *)realloc(value, len);
    memcpy(value, ptr, len);
    *(value + len) = '\0';
    pos3 = (double) atof(value);
    ptr++;

    // Ensure we get it right
    for (_ptr = ptr; *_ptr != '\n'; _ptr++, len++)
            ;

    // Store remaining atoms for future
    _ptr++;
    *fileBuffer = _ptr;

    // Now we have all we need and need to insert the numbers where they should go

    atom1->CHARGE = atom_charge;
    atom1->POS1 = pos1;
    atom1->POS2 = pos2;
    atom1->POS3 = pos3;

    // Check if next atom is worth considering
    if (*_ptr == '#')
    {
        atom1->CHECK_NEXT = 0;
    }
}


int main()
{
    // Definerer noen variabler, ikke bry deg om at de er veldig rare:
    int n_Nuclei=5000,/* redefined later */ n_Electrons; bool Test_Ongoing, Relax_Pos;
    string Basis_Set, Method; double Energy, convergance_criteria;
    vec Z = zeros(n_Nuclei); mat R = zeros(n_Nuclei, 3); // redefined later
    clock_t start = clock(); bool print_stuffies; one_atom ATOMM;
    clock_t slutt; int rank, size; bool use_angstrom; int i;
    ATOMM.CHECK_NEXT = 1; // This is used to check if there are more atoms, this way we dont have to define how meny atoms we are using in input file
    // This is done to simplify input

    // Gogo MPI
    int argc;
    char **argv;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // We want to read input from a seperate text file named INCAR, this is oldschool C (in my perspective)

    // /*
        int filefd = 0;
        struct stat file_buf;
        if ( (filefd = open("INCAR", O_RDONLY)) == -1)
        {
            perror("Failed to open file");
            exit(1);
        }
        if (fstat(filefd, &file_buf) == -1)
        {
            perror("Failed to open file");
            exit(1);
        }
        char buff[file_buf.st_size];
        read(filefd, buff, file_buf.st_size);

        // Get Basis Set
        char *retval = find_and_write(buff, "Basis_Set");
        Basis_Set = retval;

        // Get Method
        retval = find_and_write(buff, "Method");
        Method = retval;

        // Get Convergance Criteria
        retval = find_and_write(buff, "convergance_criteria");
        convergance_criteria = pow(10.0, atof(retval));

        // Get Extras
        retval = find_and_write(buff, "Relax_Pos");
        if (strncmp(retval, "true", 4) == 0)
          Relax_Pos = true;
        else
          Relax_Pos = false;

        retval = find_and_write(buff, "use_angstrom");
        if (strncmp(retval, "true", 4) == 0)
          use_angstrom = true;
        else
          use_angstrom = false;

        retval = find_and_write(buff, "print_stuffies");
        if (strncmp(retval, "true", 4) == 0)
          print_stuffies = true;
        else
          print_stuffies = false;

        // Find atoms for input here
        retval = find_start_of_atoms(buff);

        i = 0;
        while (ATOMM.CHECK_NEXT == 1)
        //for (int i = 0; i < 3; i++)
        {
            Fill_One_Atom(&retval, &ATOMM);
            Z(i) = ATOMM.CHARGE;
            R(i,0) = ATOMM.POS1;
            R(i,1) = ATOMM.POS2;
            R(i,2) = ATOMM.POS3;
            i++;
        }

        // */


        // If you dont want to use INCAR, use some of this, and comment out that code above

    // Some general input, which is just basics
    //convergance_criteria = pow(10.0, -8.0); // Criteria at which calculations broken
    Test_Ongoing = false;

/*
    // Pick basis set, can be read from INCAR
    //Basis_Set = "Thjiisen";
    //Basis_Set = "STO-3G";
    //Basis_Set = "3-21G";
    //Basis_Set = "4-31G";
    //Basis_Set = "cc-pVDZ";
    //Basis_Set = "cc-pVTZ";
    //Basis_Set = "6-311-2d2p";

    // Pick method
    //Method = "CCSD";
    //Method = "HF";
*/
    // Choose some extras, none of which works jet
    //Relax_Pos = false;
    bool Big_Tit_Run = false;
    //use_angstrom = false;
    //print_stuffies = true;

    // Input posisjon og ladninger:
    //Z(0) = 10;
/*
    Z(0) = 1;
    Z(1) = 8;
    Z(2) = 1;

    R(0,0) = 0;
    R(0,1) = 1.079252144093028;
    R(0,2) = 1.474611055780858;
    R(2,0) = 0;
    R(2,1) = 1.079252144093028;
    R(2,2) = -1.474611055780858;
    R(1,0) = 0;
    R(1,1) = 0;
    R(1,2) = 0;
*/

    // Input complete, press play.

    // Lets go
    if (Test_Ongoing == false)
    {
        Initializer Init(n_Nuclei, n_Electrons, Basis_Set, Method, Relax_Pos, R, Z, convergance_criteria, rank, size, use_angstrom);
        Energy = Init.Go(print_stuffies);
    }

    if (Test_Ongoing == true)
    {
        Hartree_Integrals HartInt;
        double x;
        double n;

        x = 57.1625;
        n = 1;

        cout << (HartInt.Boys(x, n)) << endl;

    }

    // Finalize
    slutt = clock();
/*
    // We want to write to an output file OUTCAR that will contain our output
    ofstream myfile;
    myfile.open("OUTCAR");
    myfile << "Tid brukt: " << (double)(slutt - start)/CLOCKS_PER_SEC << "s" << endl;
    myfile << "Energy: " << Energy << endl;
    myfile.close();
*/
    if (print_stuffies == true && rank == 0)
    {
        cout << "Tid brukt: " << (double)(slutt - start)/CLOCKS_PER_SEC << "s" << endl;
    }

    MPI_Finalize();
    return 0;
}

int return_charge(string Atom_Type)
{
    if (Atom_Type == "H")
    {
        return 1;
    }
    else if (Atom_Type == "He")
    {
        return 2;
    }
    else if (Atom_Type == "Li")
    {
        return 3;
    }
    else if (Atom_Type == "Be")
    {
        return 4;
    }
    else if (Atom_Type == "B")
    {
        return 5;
    }
    else if (Atom_Type == "C")
    {
        return 6;
    }
    else if (Atom_Type == "N")
    {
        return 7;
    }
    else if (Atom_Type == "O")
    {
        return 8;
    }
    else if (Atom_Type == "F")
    {
        return 9;
    }
    else if (Atom_Type == "Ne")
    {
        return 10;
    }
    else if (Atom_Type == "Na")
    {
        return 11;
    }
    else if (Atom_Type == "Mg")
    {
        return 12;
    }
    else if (Atom_Type == "Al")
    {
        return 13;
    }
    else if (Atom_Type == "Si")
    {
        return 14;
    }
    else if (Atom_Type == "P")
    {
        return 15;
    }
    else if (Atom_Type == "S")
    {
        return 16;
    }
    else if (Atom_Type == "Cl")
    {
        return 17;
    }
    else if (Atom_Type == "Ar")
    {
        return 18;
    }
    else if (Atom_Type == "K")
    {
        return 19;
    }
    else if (Atom_Type == "Ca")
    {
        return 20;
    }
    else if (Atom_Type == "Sc")
    {
        return 21;
    }
    else if (Atom_Type == "Ti")
    {
        return 22;
    }
    else if (Atom_Type == "V")
    {
        return 23;
    }

    else if (Atom_Type == "Cr")
    {
        return 24;
    }
    else if (Atom_Type == "Mn")
    {
        return 25;
    }
    else if (Atom_Type == "Fe")
    {
        return 26;
    }
    else if (Atom_Type == "Co")
    {
        return 27;
    }
    else if (Atom_Type == "Ni")
    {
        return 28;
    }
    else if (Atom_Type == "Cu")
    {
        return 29;
    }
    else if (Atom_Type == "Zn")
    {
        return 30;
    }
    else if (Atom_Type == "Ga")
    {
        return 31;
    }
    else if (Atom_Type == "Ge")
    {
        return 32;
    }

    else if (Atom_Type == "As")
    {
        return 33;
    }
    else if (Atom_Type == "Se")
    {
        return 34;
    }
    else if (Atom_Type == "Br")
    {
        return 35;
    }
    else if (Atom_Type == "Kr")
    {
        return 36;
    }
    else if (Atom_Type == "Rb")
    {
        return 37;
    }
    else if (Atom_Type == "Sr")
    {
        return 38;
    }
    else if (Atom_Type == "Y")
    {
        return 39;
    }
    else if (Atom_Type == "Zr")
    {
        return 40;
    }
    else if (Atom_Type == "Nb")
    {
        return 41;
    }

    else if (Atom_Type == "Mo")
    {
        return 42;
    }
    else if (Atom_Type == "Tc")
    {
        return 43;
    }
    else if (Atom_Type == "Ru")
    {
        return 44;
    }
    else if (Atom_Type == "Rh")
    {
        return 45;
    }
    else if (Atom_Type == "Pd")
    {
        return 46;
    }
    else if (Atom_Type == "Ag")
    {
        return 47;
    }
    else if (Atom_Type == "Cd")
    {
        return 48;
    }
    else if (Atom_Type == "In")
    {
        return 49;
    }
    else if (Atom_Type == "Sn")
    {
        return 50;
    }
    else if (Atom_Type == "Sb")
    {
        return 51;
    }
    else if (Atom_Type == "Te")
    {
        return 52;
    }
    else if (Atom_Type == "I")
    {
        return 53;
    }
    else if (Atom_Type == "Xe")
    {
        return 54;
    }
    else if (Atom_Type == "Cs")
    {
        return 55;
    }
    else if (Atom_Type == "Ba")
    {
        return 56;
    }
    else if (Atom_Type == "La")
    {
        return 57;
    }
    else if (Atom_Type == "Hf")
    {
        return 58;
    }
    else if (Atom_Type == "Ta")
    {
        return 59;
    }
    else if (Atom_Type == "W")
    {
        return 60;
    }
    else if (Atom_Type == "Re")
    {
        return 61;
    }
    else if (Atom_Type == "Os")
    {
        return 62;
    }
    else if (Atom_Type == "Ir")
    {
        return 63;
    }
    else if (Atom_Type == "Pt")
    {
        return 64;
    }
    else if (Atom_Type == "Au")
    {
        return 65;
    }
    else if (Atom_Type == "Hg")
    {
        return 66;
    }
    else if (Atom_Type == "Tl")
    {
        return 67;
    }
    else if (Atom_Type == "Pb")
    {
        return 68;
    }
    else if (Atom_Type == "Bi")
    {
        return 69;
    }
    else if (Atom_Type == "Po")
    {
        return 70;
    }
    else if (Atom_Type == "At")
    {
        return 71;
    }
    else if (Atom_Type == "Rn")
    {
        return 72;
    }
    else if (Atom_Type == "Fr")
    {
        return 73;
    }
    else if (Atom_Type == "Ra")
    {
        return 74;
    }
    else if (Atom_Type == "Ac")
    {
        return 75;
    }
    else if (Atom_Type == "Rf")
    {
        return 76;
    }
    else if (Atom_Type == "Db")
    {
        return 77;
    }
    else if (Atom_Type == "Sg")
    {
        return 78;
    }
    else if (Atom_Type == "Bh")
    {
        return 79;
    }
    else if (Atom_Type == "Hs")
    {
        return 80;
    }
    else if (Atom_Type == "Mt")
    {
        return 81;
    }
    else if (Atom_Type == "Uun")
    {
        return 82;
    }
    else if (Atom_Type == "Uuu")
    {
        return 83;
    }
    else if (Atom_Type == "Uub")
    {
        return 84;
    }
    else if (Atom_Type == "Uut")
    {
        return 85;
    }
    else if (Atom_Type == "Uuq")
    {
        return 86;
    }
    else if (Atom_Type == "Uup")
    {
        return 87;
    }
    else if (Atom_Type == "Uuh")
    {
        return 88;
    }
    else if (Atom_Type == "Uus")
    {
        return 89;
    }
    else if (Atom_Type == "Uuo")
    {
        return 90;
    }
    else if (Atom_Type == "Ce")
    {
        return 91;
    }
    else if (Atom_Type == "Pr")
    {
        return 92;
    }
    else if (Atom_Type == "Nd")
    {
        return 93;
    }
    else if (Atom_Type == "Pm")
    {
        return 94;
    }
    else if (Atom_Type == "Sm")
    {
        return 95;
    }
    else if (Atom_Type == "Eu")
    {
        return 96;
    }
    else if (Atom_Type == "Gd")
    {
        return 97;
    }
    else if (Atom_Type == "Tb")
    {
        return 98;
    }
    else if (Atom_Type == "Dy")
    {
        return 99;
    }
    else if (Atom_Type == "Ho")
    {
        return 100;
    }
    else if (Atom_Type == "Er")
    {
        return 101;
    }
    else if (Atom_Type == "Tm")
    {
        return 102;
    }
    else if (Atom_Type == "Yb")
    {
        return 103;
    }
    else if (Atom_Type == "Lu")
    {
        return 104;
    }
    else if (Atom_Type == "Th")
    {
        return 105;
    }
    else if (Atom_Type == "Pa")
    {
        return 106;
    }
    else if (Atom_Type == "U")
    {
        return 107;
    }
    else if (Atom_Type == "Np")
    {
        return 108;
    }
    else if (Atom_Type == "Pu")
    {
        return 109;
    }
    else if (Atom_Type == "Am")
    {
        return 110;
    }
    else if (Atom_Type == "Cm")
    {
        return 111;
    }
    else if (Atom_Type == "Bk")
    {
        return 112;
    }
    else if (Atom_Type == "Cf")
    {
        return 113;
    }
    else if (Atom_Type == "Es")
    {
        return 114;
    }
    else if (Atom_Type == "Fm")
    {
        return 115;
    }
    else if (Atom_Type == "Md")
    {
        return 116;
    }
    else if (Atom_Type == "No")
    {
        return 117;
    }
    else if (Atom_Type == "Lr")
    {
        return 118;
    }

}
