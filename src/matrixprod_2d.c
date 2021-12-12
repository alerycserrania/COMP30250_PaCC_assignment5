#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>


void print_blocked_matrix(int N, int BS, double* M)
{
    for (int bsi = 0; bsi < N; bsi+=BS)
    {
        for (int bi = 0; bi < BS; bi++)
        {
            for (int bsj = 0; bsj < N; bsj+=BS)
            {
                for (int bj = 0; bj < BS; bj++)
                {
                    printf("%f\t ", M[bsi*BS + bi*BS + bsj*N + bj]);
                }
            }
            printf("\n");
        }
    }
}


int main(int argc, char *argv[])
{
    int m_size, m_part, nb_procs, rank;
    double *a, *b, *c;
    double *all_a, *all_b, *all_c;
    double *arows, *bcols;
    int row, col;
    int i, j, k, ind1, ind2;
    double start, sum;

    MPI_Comm rowcomm, colcomm;

    char PRINT_RESULT, PRINT_TIME;

    if(argc != 4)
    {
       printf("Please, use: %s N PRINT_RESULT PRINT_TIME:\n", argv[0]);
       printf("\t- N: matrix size\n");
       printf("\t- PRINT_RESULT (y/n): print result to stdout\n");
       printf("\t- PRINT_TIME (y/n): print time elasped to stdout\n");
       exit(EXIT_FAILURE);
    }

    m_size = atoi(argv[1]);
    PRINT_RESULT = argv[2][0];
    PRINT_TIME = argv[3][0];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    m_part = m_size/(int)sqrt(nb_procs);

    a = malloc(m_part * m_part * sizeof(double));
    b = malloc(m_part * m_part * sizeof(double));
    c = malloc(m_part * m_part * sizeof(double));
    arows = malloc(m_size * m_part * sizeof(double));
    bcols = malloc(m_size * m_part * sizeof(double));
    row = rank / (m_size/m_part);
    col = rank % (m_size/m_part);

    for (i = 0; i < m_part * m_part; i++)
    {
        a[i] = 1. + i;
        b[i] = 1. + i;
        c[i] = 0.;
    }

    MPI_Comm_split(MPI_COMM_WORLD, row, rank, &rowcomm);
    MPI_Comm_split(MPI_COMM_WORLD, col, rank, &colcomm);

    if (rank == 0)
    {
        start = MPI_Wtime();
    }

    MPI_Allgather(a, m_part * m_part, MPI_DOUBLE, arows, m_part * m_part, MPI_DOUBLE, rowcomm);
    MPI_Allgather(b, m_part * m_part, MPI_DOUBLE, bcols, m_part * m_part, MPI_DOUBLE, colcomm);

    for (i = 0; i < m_part; i++)
    {
        for (j = 0; j < m_part; j++)
        {
            sum = 0.0;
            ind1 = 0;
            ind2 = 0;
            for (k = 0; k < m_size; k++)
            {
                sum += arows[ind1+ind2+i*m_part] * bcols[m_part*k+j];
                if (ind1 < m_part-1) {
                    ind1 += 1;
                } else {
                    ind1 = 0;
                    ind2 += m_part*m_part;
                }
            }
            c[i*m_part+j] = sum;
        }
    }

    if (rank == 0 && PRINT_TIME == 'y')
    {
        printf("elapsed: %fs\n", MPI_Wtime() - start);
    }

    if (PRINT_RESULT == 'y')
    {
        if (rank == 0) 
        {
            all_a = malloc(m_size * m_size * sizeof(double));
            all_b = malloc(m_size * m_size * sizeof(double));
            all_c = malloc(m_size * m_size * sizeof(double));
        }

        MPI_Gather(a, m_part * m_part, MPI_DOUBLE, all_a, m_part * m_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(b, m_part * m_part, MPI_DOUBLE, all_b, m_part * m_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(c, m_part * m_part, MPI_DOUBLE, all_c, m_part * m_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (rank == 0)
        {
            printf("A =\n");
            print_blocked_matrix(m_size, m_part, all_a);
            printf("B =\n");
            print_blocked_matrix(m_size, m_part, all_b);
            printf("C =\n");
            print_blocked_matrix(m_size, m_part, all_c);
        }
    }



    MPI_Comm_free(&rowcomm);
    MPI_Comm_free(&colcomm);

    free(a);
    free(b);
    free(c);
    free(arows);
    free(bcols);

    MPI_Finalize();

    return 0;
}

