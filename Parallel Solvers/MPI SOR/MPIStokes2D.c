/*
    NOTE: Parts of the solver are withheld but available on request.  Structure, 
    documentation, and the full algorithm description are intact.  This file
    will not compile as posted.
*/

/*
    Solve the Stokes equations using a MAC grid that is approximately N by N.
    The u array should be N by (N - 1), the v array should be (N - 1) by N,
    and the p array should be (N - 1) by (N - 1).  Experimentation will be used
    to find the omega value that requires the fewest iteration to reach a 
    tolerance of 10^-9.  This will be done using MPI with the foundation modified
    from the OpenMP version.

    Stencils are being used for the partial derivatives in the x and y direction 
    but the residuals for the three arrays around found and updated in place.  
    The maximum residual over each sweep is used to determine the convergence.
    Careful attention is given to satisfy the boundary conditions at each step.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

/*
int main(int argc, char* argv[])

    This is the main function that uses SOR to solve 2D stokes flow in a region 
    with no starting flow, boundaries on the top and bottom, and a pressure of 1 
    on the right.

    Inputs:
        int argc:   This is the number of parameters passed to the program, this 
        function can take 6 values.

        char* argv[]:   The array of parameters passed to the program, the 
            description of each is below.
            
            argv[1]:    int N = grid size
            argv[2]:    double mu = viscosity
            argv[3]:    double P = pressure drop
            argv[4]:    double omega = relaxation parameter
            argv[5]:    double tolerance = error tolerance
            argv[6]:    int K = maximum number of iterations

    Outputs:
        int:    0 if successful and nonzero if it was unsuccessful
        StokesX.out:    A binary output file that stores the results of the 
                        SOR solver, X can be U, V, or P
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
    // to have memory allocated.  The problem has 6 parameters or arguments
    //      N = grid size
    //      mu = viscosity
    //      P = pressure drop
    //      omega = relaxation parameter
    //      tolerance = error tolerance
    //      K = maximum number of iterations
    int N;
    double mu;
    double P;
    double omega;
    double tolerance;
    int K;

    // Have the first process determine the values to be broadcast
    if (rank == 0)
    {
        // Read the values from the argument list if given, not ideal will do better
        // version for next set of projects
        if (argc == 7)
        {
            N = atoi(argv[1]);
            mu = atof(argv[2]);
            P = atof(argv[3]);
            omega = atof(argv[4]);
            tolerance = atof(argv[5]);
            K = atoi(argv[6]);
        }

        // Otherwise use default values
        else
        {
            N = 128;
            mu = 1.0;
            P = 1.0;
            omega = 0.47;
            tolerance = 1e-9;
            K = 100000;
        }
    }

    // Broadcast the values
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&P, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*
        For this problem I am going to treat the three arrays as 1 N x N array in terms 
        of how I am going to split the arrays between the cores.  Depending on the 
        grid size and number of cores used I am going to split the rows as evenly as 
        possible, with the earlier ranks receiving more rows if there is a remainder.
        For example if the first rank is assigned rows 0 through 5 then it is going to 
        work with this subset of rows for the u, v, and p arrays.  While this is how 
        the rows are allocated it is important to remember that the rows themselves 
        vary between arrays and that the u and p arrays have N - 1 rows.

        The other important factor to consider is that there needs to be additional
        values stored to account for the stencil needed values in the rows above and 
        below.  For the first rank this means that there will be one extra row below,
        for the last rank this means there will be one extra row above, and for the
        rest of the ranks it means that there will be an extra row above and below the
        rows being updated.

        Finally, there will be two points where the data is shared between the cores 
        at each iteration.  The first is after the red points have been updated in
        each array and the second is after the black points have been updated in each
        array.  At both of these points, the shared rows will be shared between the
        necessary cores.  My first implementation will involve a transfer for each
        array but it might be possible to transfer the data as a struct.  After the 
        solution has converged all of the portions will be gathered and saved.  

        The last part to consider with the SOR solver is how the max residual will
        be shared between cores.  This should be easier than my OpenMP implementation
        since I only need to share it once, after the black points were updated.
    */

    /*
        The number of rows in the processes solution arrays depends on the rank of the
        CPU, the first and last rank have one fewer row since there is only one padding
        row needed and the last rank also have 1 fewer rows for the u and p arrays due
        to the structure of the MAC grid.  
        
        The total number of rows, number of solution rows, and true first row in the 
        larger array will be stored for each process for rebuilding the solution.  In 
        addition to these value each process will store the size, number of columns, and
        number of rows for each of the sub u, v, and p arrays.
    */

    // REDACTED: array initialization and communication values are available on request

    /*
        Non blocking communication is going to be used to transfer data between processes
        so the status of these transfers needs to be tracked.  To do this there will be 
        pairs of vectors and sets of MPI_Request objects for each array.  The first
        vector will track if the data is ready as it cannot be used while a transfer is
        in progress, and the second will track if the edge rows have been updated to
        ensure it is not done multiple times.  There will be an MPI_Request for the top
        and bottom send and receives.  The structure for the data transfer vector will
        take the following form with x = {u, v, p},
            xReady[0] = top send status
            xReady[1] = top request status
            xReady[2] = bottom send status
            xReady[3] = bottom request status
        The update status vector will have this form,
            xDone[0] = lowest index row updated
            xDone[1] = highest index row updated
        For both of the vectors, a 1 will represent the transfer being finished or if
        the updates were complete.  Since the value are set to 1, the values will be
        reset at the start of each step and will allow for easier handling of the edge
        cases.  Also, to make conditional logic easier, all of the edge cases will be
        handled at the same time so the done vectors will always reset so both values
        are 0.  ** The top send is at the top of the array to the (rank - 1) process
        and the bottom send is at the bottom of the array to the (rank + 1) process. **
    */

    // // U array
    MPI_Request uSendTopRequest;
    MPI_Request uSendBottomRequest;
    MPI_Request uRecieveTopRequest;
    MPI_Request uRecieveBottomRequest;

    int uReady[4] = {1, 1, 1, 1};

    int uDone[2] = {1, 1};

    // // V array
    MPI_Request vSendTopRequest;
    MPI_Request vSendBottomRequest;
    MPI_Request vRecieveTopRequest;
    MPI_Request vRecieveBottomRequest;

    int vReady[4] = {1, 1, 1, 1};

    int vDone[2] = {1, 1};

    // // P array
    MPI_Request pSendTopRequest;
    MPI_Request pSendBottomRequest;
    MPI_Request pRecieveTopRequest;
    MPI_Request pRecieveBottomRequest;

    int pReady[4] = {1, 1, 1, 1};

    int pDone[2] = {1, 1};

    // This is going to track the number of arrays finished and which ones are 
    // finished.  The count will range from 0 to 3 and this will be used to
    // in the conditional in the while loop to check flags
    //      arrayStatus[0] = u array status
    //      arrayStatus[1] = v array status 
    //      arrayStatus[2] = p array status
    int arrayStatus[3] = {1, 1, 1};
    int finishedArrays = 3;

    // Pre calculate some of the coefficients, dx = dy = 1 / (N - 1) so dx / dy
    // = dy / dx = 1 which will be assumed in the SOR loop.
    double dx_dy = 1.0 / (N - 1);

    /*
        For the SOR loop, all of the processes need access to the maximum residual
        and iteration count to know when to exit the while loop.  This is going to be
        controlled by the rank 0 process with communication happening at the end of the
        updates.  A do while loop is used since the value of the maximum residual is not
        known before and this guarantees at least 1 step of the algorithm.

        The structure of the red black SOR algorithm used is:
            1) Update all of the red points
            2) Transfer the updated values
            3) Update the middle black points (rows that don't rely on transferred data)
            4) Reset flags + variables that track transfer status
            5) Use flags to wait for transfers to finish and update the rest of the black
                points as data completes transfering
            6) Transfer the updated values
            7) Wait for transfer to complete (MPI_Wait)
            8) Find the max residual with MPI_Allreduce
    */
    double residual = 0.0;
    double localMaxResidual = 0.0;
    double maxResidual = 0.0;
    int iterationNumber = 0;

    // REDACTED: SOR do while loop is available on request

    // Have all processes end the timer and then determine the max time
    double elapsedTime = MPI_Wtime() - startTime;
    double maxElapsedTime;
    MPI_Allreduce(&elapsedTime, &maxElapsedTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /*
        To collect the solution using All_Gatherv, the first ranks needs to define three
        pairs of vectors to store the number of values from each process and the correct
        displacement for those values.  Afterwards All_Gatherv can be called once per
        array to create the final solution.
    */

    // REDACTED: final solution code is available on request

    // Rank 0 writes the solution to memory and prints the timing results
    if (rank == 0)
    {
        // Open the files that will store the data
        FILE* fileid_u = fopen("StokesU.out", "w");
        FILE* fileid_v = fopen("StokesV.out", "w");
        FILE* fileid_p = fopen("StokesP.out", "w");

        // Save the data
        fwrite(finalU, sizeof(double), (N * (N - 1)), fileid_u);
        fwrite(finalV, sizeof(double), (N * (N - 1)), fileid_v);
        fwrite(finalP, sizeof(double), ((N - 1) * (N - 1)), fileid_p);

        // Close the files used
        fclose(fileid_u);
        fclose(fileid_v);
        fclose(fileid_p);

        free(finalU);
        free(finalV);
        free(finalP);
        free(uCounts);
        free(uDisplacements);
        free(vCounts);
        free(vDisplacements);
        free(pCounts);
        free(pDisplacements);

        // Print parameters and results
        double timePerIteration = maxElapsedTime / iterationNumber;
        printf("Grid Length: %d\n", N);
        printf("Omega Value: %f\n", omega);
        printf("Number of Cores: %d\n", size);
        printf("Number of Iterations: %d\n", iterationNumber);
        printf("Function took %f seconds to execute\n", maxElapsedTime);
        printf("Average Time per Iteration: %lf\n", timePerIteration);
        printf("Max Ending Residual: %.17f\n", maxResidual);
    }

    free(u);
    free(v);
    free(p);
    
    MPI_Finalize();

    return 0;
}