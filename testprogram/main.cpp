#include <stdio.h>
#include "string.h"
#include <stdlib.h>
#include <math.h>
#include "../../home/ole/Desktop/include/armadillo"
#include <iostream>

using namespace arma;
using namespace std;

/*
    This program is all about parallel calculation of matrix multiplications

    We will use both MPI and OpenMP, but seperatly. Not in combination




*/


// Pre defined functions
void read_matrix_binaryformat (char* filename, double*** matrix, int* num_rows, int* num_cols);
void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols);
void allocate_matrix(double*** matrix,int* num_rows,int* num_cols);


int main(void)
{
    // Declarations
    double **a_matrix;
    double **b_matrix;
    double **c_matrix;
    int a_x, a_y;
    int b_x, b_y;
    int i, j, k, l;

    mat A = zeros(100, 50);
    mat B = zeros(50, 100);


    // Initialization
    read_matrix_binaryformat("small_matrix_a.bin", &a_matrix, &a_x, &a_y);
    read_matrix_binaryformat("small_matrix_b.bin", &b_matrix, &b_x, &b_y);
    allocate_matrix(&c_matrix, &a_x, &b_y);


    for (i = 0; i < a_x; i++)
    {
        for (j=0; j< a_y; j++)
        {
            A(i,j) = a_matrix[i][j];
            B(j,i) = b_matrix[j][i];
        }
    }

    // Initialize c_matrix to 0
    for (i=0; i<a_x; i++)
    {
        for (j=0; j<b_y; j++)
        {
            c_matrix[i][j] = 0;
        }
    }

    // Calculate matrix c
    for (i=0; i<a_x; i++)
    {
        for (j=0; j<b_y; j++)
        {
            for (k=0; k<a_y; k++)
            {
                c_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
            }
        }
    }


    mat C = A * B;

    // Print matrix a
    for (i=0; i<a_x; i++)
    {
        for (l=0; l<b_y; l++)
        {
            printf("%f ", c_matrix[i][l]);
            cout << C(i,l) << endl;
        }
    }

    printf("%f", c_matrix[0][0]);


    // Write matrix
    //write_matrix_binaryformat(output_matrix_name, matrix, num_x, num_y);
    // Finalization
    printf("\n Hello World!\n");
    return 0;
}


void read_matrix_binaryformat(char* filename,double*** matrix,int* num_rows, int* num_cols)
{

  size_t result;
  FILE* fp = fopen(filename,"rb");
  result = fread(num_rows,sizeof(int),1,fp);
  result = fread(num_cols,sizeof(int),1,fp);
  printf("Loading matrix %d x %d \n",*num_rows,*num_cols);
  allocate_matrix(matrix,num_rows,num_cols);
  result = fread((*matrix)[0],sizeof(double),(*num_rows)*(*num_cols),fp);
  fclose(fp);
}


void allocate_matrix(double*** matrix,int* num_rows,int* num_cols)
{

    int i;

    *matrix = (double**) malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*) malloc((*num_rows)*(*num_cols)*sizeof(double));

    for (i = 1 ; i < (*num_rows); i++)
    {
        (*matrix)[i] = (*matrix)[i-1] + (*num_cols);
    }
}



void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&num_rows, sizeof(int), 1, fp);
    fwrite (&num_cols, sizeof(int), 1, fp);
    fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
    fclose (fp);
}
