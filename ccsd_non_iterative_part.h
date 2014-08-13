#ifndef CCSD_NON_ITERATIVE_PART_H
#define CCSD_NON_ITERATIVE_PART_H

#include "../../home/ole/Desktop/include/armadillo"
#include <iostream>
#include "mpi.h"

using namespace std;
using namespace arma;

class ccsd_non_iterative_part
{
public:
    ccsd_non_iterative_part(int unc, int elec, int ran, int siz);
    int rank, size, n_Electrons, unocc_orb;
    vec WORK_EACH_NODE;
    vec WORK_EACH_NODE_Part1;
    vec WORK_EACH_NODE_Part2;
    int** Global_Displacement_1;
    int** Global_Displacement_2;
    int** Global_Worksize_1;
    int** Global_Worksize_2;
    int max_terms;

    int Local_Displacement_1;
    int Local_Displacement_2;
    int Local_Displacement_3;

    int Calc_sum_a_n(int a);

    int ret_max_terms();

    void Map_T2_For_MPI();
    vec return_Work_T2();

    void Map_Part1_For_MPI();
    vec return_Work_P1();
    vec return_Work_P2();
    int** ret_Global_Disp1();
    int** ret_Global_Disp2();
    int** ret_Global_Work1();
    int** ret_Global_Work2();

    mat Where_To_Start_Part1;
    mat Where_To_Start_Part2;
    mat Where_To_start_T2;

    void Figure_Out_Where_To_Start();
    mat Return_Start_Part2_Pos();
    mat Return_Start_T2_Pos();
};

#endif // CCSD_NON_ITERATIVE_PART_H
