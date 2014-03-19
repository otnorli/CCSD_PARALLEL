#include "matrix_size_setter.h"

matrix_size_setter::matrix_size_setter(vec zZz, string sys_string, int n_N)
{
    System_String = sys_string;
    Z = zZz;
    n_Nuclei = n_N;
}

int matrix_size_setter::Set_Matrix_Size()
{
    uint i;
    value = 0;

    for (i=0; i<n_Nuclei; i++)
    {
        if (System_String == "3-21G")
        {
            if (Z(i) == 1)
            {
                value += 2;
            }

            if (Z(i) == 4)
            {
                value += 5;
            }

            // Flere atomer
        }

        if (System_String == "4-31G")
        {
            if (Z(i) == 8)
            {
                value += 5;
            }
        }

        if (System_String == "6-311G")
        {
            if (Z(i) == 1)
            {
                value += 3;
            }
        }

        if (System_String == "Thjiisen")
        {
            if (Z(i) == 1)
            {
                value += 4;
            }
        }
    }

    return value;
}
