/*
    NOTE: Parts of the solver are withheld but available on request.  Structure, 
    documentation, and the full algorithm description are intact.  This file
    will not compile as posted.
*/

/* 
    This is the MPI implementation of the linearized Euler equations where
    each process has the full array.  An attempt was made to use scatter 
    calls but this was unsuccessful
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

extern void dgetrf_(int*, int*, double*, int*, int*, int*);
extern void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);

/*
    This is the main function that uses ADI to solve the linearized Euler 
    equations.

    Inputs:
        int argc:   This is the number of parameters passed to the program, this function can take
                    2 values.

        char* argv[]:   The array of parameters passed to the program, the description of each is 
                        below.
            argv[1]:    int N = grid size
            argv[2]:    int M = number of time steps

    Outputs:
        int:    0 if successful and nonzero if it was unsuccessful
        StokesX.out:    A binary output file that stores the results of the SOR solver,
                        X can be U, V, or P
*/
int main(int argc, char* argv[])
{
    // Initialize the MPI communication
    MPI_Init(&argc, &argv);

    // Handle error checking using the header file, run with
    //      mpicc -c mpiprog.c -DDO_ERROR_CHECKING
    #ifdef DO_ERROR_CHECKING
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    #endif

    // Determine the rank of the current process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine the total number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Have all processes call the timing function
    double precision = MPI_Wtick();
    double startTime = MPI_Wtime();

    // The first process will broadcast the parameters but each process needs
    // to have memory allocated for the parameters and other constants
    //      N = grid size
    //      M = time steps
    int N;
    int M;
    double gamma;
    double dt;
    double dN;

    // Have the first process determine the values to be broadcast
    if (rank == 0)
    {
        // Read the values from the argument list if given, not ideal will do better
        // version for next set of projects
        if (argc == 3)
        {
            N = atoi(argv[1]);
            M = atoi(argv[2]);
        }

        // Otherwise use default values
        else
        {
            N = 128;
            M = 2000;
        }
    }

    // Broadcast the values
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Set the other default values
    gamma = 1.4;
    dt = 2.0 / M;
    dN = 2.0 / N;

    /*
        For my implementation of ADI each grid will store the full array for each
        of the variables.  As a result Allgatherv calls will be used since I am not
        requiring that the grid size be divisible bu the number of processes used.
        For these data transfers to work vectors need to be created that store the
        number of rows each process is working with, the starting offset, and the
        starting row.
    */

    // REDACTED: offset and row count code available on request

    /*
        Each process is going to get the full set of rows and populate them using
        the following initial conditions
            row(x, y, 0) = 2/gamma * e^{-100(x^2 + y^2)}
            u(x, y, 0) = 0
            v(x, y, 0) = 0
            p(x, y, 0) = 2e^{-100(x^2 + y^2)}
        This will not be done in parallel, each process will populate the initial
        conditions.  Two sets of arrays will be needed since updates are not in 
        place, one for the 'old' values and one for the 'new' values
    */
    double (*rowOld)[N] = (double(*)[N])malloc((N) * sizeof(*rowOld));
    double (*uOld)[N] = (double(*)[N])malloc((N) * sizeof(*uOld));
    double (*vOld)[N] = (double(*)[N])malloc((N) * sizeof(*vOld));
    double (*pOld)[N] = (double(*)[N])malloc((N) * sizeof(*pOld));

    double (*rowNew)[N] = (double(*)[N])malloc((N) * sizeof(*rowNew));
    double (*uNew)[N] = (double(*)[N])malloc((N) * sizeof(*uNew));
    double (*vNew)[N] = (double(*)[N])malloc((N) * sizeof(*vNew));
    double (*pNew)[N] = (double(*)[N])malloc((N) * sizeof(*pNew));

    // REDACTED: initialization code available on request

    /*
        For the implicit steps we are using 2 matrices for each step of the solver.
        These do not change so I am going initialize them and prefactor them to save
        time during the ADI loop.  These matrices are 4N by 4N and have the following
        block structure (each block is N by N):
                          | I   0   D   0 |
              implicitY = | 0   I   0   0 |
                          | 0   0   I   D |
                          | 0   0  gD   I |

                          | I   D   0   0 |
              implicitX = | 0   I   0   D |
                          | 0   0   I   0 |
                          | 0  gD   0   I |
        Where I is the identity matrix, D which has -dt/(4 * dx) on the sub and super
        diagonal, and gD being gamma * D.  For D and gD there is wrap around with the
        sub and super values that push them to the opposite edge when they go out of
        bounds.  To keep track of the current block during initialization I am going 
        to keep track of the current row block and column block with these values 
        being 1 indexed.
    */
    int rowBlock = 0;
    int colBlock = 0;

    // REDACTED: solution initialization code available on request

    /*
        Each process will factor the implicit matrices using LAPACK, scaLAPACK is
        not being used since the parallelism comes from splitting the updates across
        multiple processes.  As a result the next step is to transpose the matrices
        into column major form and then use dgetrf.
    */

   // REDACTED: factorization code is available on request
    
    // Prepare the files for saving the results
    FILE* fileid_row = NULL;
    FILE* fileid_u = NULL;
    FILE* fileid_v = NULL;
    FILE* fileid_p = NULL;
    
    // Have rank 0 save the initial results
    if (rank == 0)
    {
        // Open the files that will store the data
        fileid_row = fopen("EulerR.out", "a");
        fileid_u = fopen("EulerU.out", "a");
        fileid_v = fopen("EulerV.out", "a");
        fileid_p = fopen("EulerP.out", "a");

        // Save the initial data
        fwrite(rowOld, sizeof(double), (N * N), fileid_row);
        fwrite(uOld, sizeof(double), (N * N), fileid_u);
        fwrite(vOld, sizeof(double), (N * N), fileid_v);
        fwrite(pOld, sizeof(double), (N * N), fileid_p);
    }

    /*
        The last step before the algorithm is to allocate memory for the RHS of
        the implicit steps.  This stores information about all 4 arrays so there 
        will be $N rows and localRows number of columns as each process might have
        a different number of assigned rows.  For this implementation, the RHS
        will be pre transposed into column major form for LAPACK.
            rightHandSideT = fourN columns by localRows rows 
    */
    nrhs = localRows;
    double (*rightHandSideT)[fourN] = (double(*)[fourN])malloc(localRows * sizeof(*rightHandSideT));

    /*
        The general form of the ADI algorithm is:
            1) Explicit X update
            2) Implicit Y update
            3) Explicit Y update
            4) Implicit X update
        During the explicit steps careful attention is given to the boundary conditions
        to ensure they are satisfied.  In addition the location of the transpose call is
        important since the C program is row major order but all calls using LAPACK are
        column major.  Lastly, the MPI gather calls are only after the explicit updates. 
    */

    // REDACTED: ADI algorithm is available on request

    // Have all processes end the timer and then determine the max time
    double elapsedTime = MPI_Wtime() - startTime;
    double maxElapsedTime;
    MPI_Allreduce(&elapsedTime, &maxElapsedTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Have the first rank close the files, free the memory and save the results
    if (rank == 0)
    {
        // Close the files used
        fclose(fileid_row);
        fclose(fileid_u);
        fclose(fileid_v);
        fclose(fileid_p);

        // Print parameters and results
        double timePerIteration = maxElapsedTime / M;
        printf("Grid Length: %d\n", N);
        printf("Time Steps: %d\n", M);
        printf("Number of Cores: %d\n", size);
        printf("Function took %f seconds to execute\n", maxElapsedTime);
        printf("Average Time per Iteration: %lf\n", timePerIteration);
    }

    // Free the arrays
    free(rankLocalRows);
    free(counts);
    free(displacements);
    free(rowOld);
    free(uOld);
    free(vOld);
    free(pOld);
    free(rowNew);
    free(uNew);
    free(vNew);
    free(pNew);
    free(implicitX);
    free(implicitY);
    free(implicitXT);
    free(implicitYT);
    free(ipivX);
    free(ipivY);
    free(rightHandSideT);
    
    MPI_Finalize();

    return 0;
}