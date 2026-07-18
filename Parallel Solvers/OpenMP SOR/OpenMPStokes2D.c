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
    tolerance of 10^-9.  This will be done using OpenMP.

    Stencils are being used for the partial derivatives in the x and y direction 
    but the residuals for the three arrays around found and updated in place.  
    The maximum residual over each sweep is used to determine the convergence.
    Careful attention is given to satisfy the boundary conditions at each step.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <float.h>

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
    // This program is going to take 6 arguments 
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

    // Read the values from the argument list if given
    if (argc == 7)
    {
        N = atoi(argv[2]);
        mu = atof(argv[3]);
        P = atof(argv[4]);
        omega = atof(argv[5]);
        tolerance = atof(argv[6]);
        K = atoi(argv[7]);
    }

    // Else set the default values
    else
    {
        N = 128;
        mu = 1.0;
        P = 1.0;
        omega = 0.4;
        tolerance = 1e-9;
        K = 100000;
    }

    // Define the number of elements in each array
    int speedSize = N * (N - 1);
    int pSize = (N - 1) * (N - 1);

    // Open the files that will store the data
    FILE* fileid_u = fopen("StokesU.out", "w");
    FILE* fileid_v = fopen("StokesV.out", "w");
    FILE* fileid_p = fopen("StokesP.out", "w");

    /*
        This solver needs 3 arrays as described above, u for the x-velocity, v
        for the y-velocity, and p for pressure.  These arrays are going to use the
        MAC structure and will use and N by N array as their base, leading to 
        each array having the following shape,
            u = N by (N - 1)
            v = (N - 1) by N
            p = (N - 1) by (N - 1)
    */
    double (*u)[N - 1] = (double(*)[N - 1])malloc(N * sizeof(*u));
    double (*v)[N] = (double(*)[N])malloc((N - 1) * sizeof(*v));
    double (*p)[N - 1] = (double(*)[N - 1])malloc((N - 1) * sizeof(*p));

    /*
        The first step is to initialize the values of the arrays to 0 which will be
        dont serially.  This section will also define the iterators that will be used
        throughout the solver, j will iterate through the x values and represent the 
        column number while k will iterate through the y values and represent the row 
        number.  Due to the row major ordering in C, k will be in the outer loop and
        j will be in the inner loop.  
    */
    int i;      // general indexer
    int j;      // j = column number        u[k][j]
    int k;      // k = row number
    int subN = N - 1;

    // Handle the majority of the arrays (N - 1) x (N - 1)
    for (k = 0; k < subN; ++k)
    {
        for (j = 0; k < subN; ++j)
        {
            u[k][j] = 0;
            v[k][j] = 0;
            p[k][j] = 0;
        }

        // Here is a good spot to set the value of the last column in the u array since
        // the memory access should be contiguous
        u[k][subN] = 0;
    }

    // A new loop will handle the last row of v to avoid memory misses in the above loop
    for (j = 0; j < subN; ++j)
    {
        v[subN][j] = 0;
    }

    /*
        One way to optimize the solver is to pre calculate the coefficients used in
        the residual calculations and then pass them to the OpenMP loop.  Also, for
        the geometry used with this problem, dx = dy = 1 / (N - 1) so dx / dy = 
        dy / dx = 1 which will be assumed in the SOR loop.
    */
    double dx_dy = 1 / (N - 1);
    double dx = 1 / (N - 1);
    double dx_2 = dx / 2;
    double dy = 1 / (N - 1);
    double dy_2 = dy / 2;
    double dxdy = dx * dy;

    /*
        Using red black ordering will allow us to utilize the parallelization that OpenMP
        gives.  The points will be split into two groups in which the updates are 
        independent and can be given to different processes.  The red points will be 
        updated first with the black point afterwards, with how C stores arrays in 
        memory this means that the red points will have an even index and the black
        points will have an odd index.  
        
        Splitting the points up this way means that the edge cases need to be handled 
        carefully leading to conditional logic for each case.  Another factor with how
        OpenMP assigns tasks is that the arrays have different sizes which is resolved
        by splitting the colored updates into three loops, one for u, v, and p.  This
        leads to the following steps:
            1) Update red points in the u array
            2) Update red points in the v array
            3) Update red points in the p array
            4) Update black points in the u array
            5) Update black points in the v array
            6) Update black points in the p array

        Since parallel structres are expensive to initialize it will be done once and 
        consist of the omp parallel block with a do while loop to check the convergence
        conditions at the end of each iteration.  Within the do while loop there are 
        parallel for loops with static scheduling and reductions that calculate the max 
        residual at the end of each step described above which tracks the largest residual 
        throughout the step.  Lastly these are the variables used in the omp chunk,
        SHARED
            mu                  viscosity
            P                   pressure difference
            omega               relaxation parameter     
            u                   x-speed array
            v                   y-speed array
            p                   pressure values
            N                   grid size
            dx_dy               dx = dy = 1 / (N - 1)
            speedSize           N * (N - 1)
            pSize               (N - 1) * (N - 1)
            numIterations       how many iterations the solver needed to converge
            maxIter...Residual  maximum residual for the entire step
    
        PRIVATE
            i                   general iterator
            j                   column number
            k                   row number
            maxResidual         max residual from that thread
            residual            residual at the current point
    */

    // Values needed for the solver
    int num_iterations = 0;
    double maxIterationResidual = -DBL_MIN;
    double maxResidual = -DBL_MIN;

    // Start the timer
    double startTime = omp_get_wtime();

    // REDACTED: OpenMP SOR loop is available on request

    // End the timer
    double endTime = omp_get_wtime();
    double cpu_time_used = endTime - startTime;
    double timePerIteration = cpu_time_used / numIterations;

    // Save the data
    fwrite(u, sizeof(double), speedSize, fileid_u);
    fwrite(v, sizeof(double), speedSize, fileid_v);
    fwrite(p, sizeof(double), pSize, fileid_p);

    // Close the files used
    fclose(fileid_u);
    fclose(fileid_v);
    fclose(fileid_p);

    // Print parameters and results
    printf("Grid Length: %d\n", N);
    printf("Number of Iterations: %d\n", numIterations);
    printf("Function took %f seconds to execute\n", cpu_time_used);
    printf("Average Time per Iteration: %f\n", timePerIteration);
    printf("Max Ending Residual: %.17f\n", maxIterationResidual);

    // Free the three arays initialized
    free(u);
    free(v);
    free(p);

    return 0;
}