/*
    NOTE: The solver body and kernels are withheld but available on request.  
    Structure, documentation, and the full algorithm description are intact.  This
    file will not compile as posted.
*/

/* 
    This is the CUDA implementation of the linearized Euler equations and
    modified from my MPI implementation.  The following kernels will be used:
        1) Data initialization
        2) Array transpose
        3) Explicit x update
        4) Explicit y update
        5) Create right hand side
        6) Disassemble right hand side
        7) Initialize implicit update arrays (?)
    Since ADI is being used and the grid size, N, is below 1024 elements, 
    each block is going to represent 1 row.

    One speedup tip given to me but I was not able to implement was treating
    the large grid as a 4x4 grid using the block structure of the matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <cuda.h>
#include <cublas_v2.h>   
#include <cusolverDn.h>

// // These are constants used in the solver
// How big should the tiles be in the transpose kernel, 16 is chosen
// instead of 32 since there will be cases when the full grid is 16 by
// 16 and it will make it more robust.  In addition this will be the 
// number of thread rows per block for this kernel
const int tile_dim = 16;

// 
__constant__ int dev_N;
__constant__ double dev_gamma;
__constant__ double dev_twoGamma;
__constant__ double dev_dN;
__constant__ double dev_dt_4_dN;

/* 
    This is a kernel to initialize the data, each row will be represented
    by one block and get the following initial conditions
        row(x, y, 0) = 2/gamma * e^{-100(x^2 + y^2)}
        u(x, y, 0) = 0
        v(x, y, 0) = 0
        p(x, y, 0) = 2e^{-100(x^2 + y^2)}

    Inputs:
        N/A

    Inputs:
        double* row:  device pointer to initialized row values
        double* u:  device pointer to initialized u values 
        double* v:  device pointer to initialized v values
        double* p:  device pointer to initialized p values
*/
__global__ void initialize(double* row, double* u, double* v, double* p)
{
    // REDACTED: kernel code is available on request
}

/* 
    This is a kernel to transpose an array in an efficient matter to avoid
    cache misses.  The tile idea was inspired from this Nvidia article,
    https://developer.nvidia.com/blog/efficient-matrix-transpose-cuda-cc/.  
    Instead of assigning the blocks to full rows like with the other kernels,
    this will split it up into square blocks to make the transpose more efficient
    in terms of memory accesses.

    Inputs:
        double* dev_source:  device pointer to source array
        int* N:  size of array

    Outputs:
        double* dev_transpose:  device pointer to transposed array
*/
__global__ void transpose(const double* dev_source, const int N, double* dev_transpose)
{
    // REDACTED: kernel code is available on request
}

/*
    This is a kernel for the explicit update in x which takes the 
    following form:
        row^(n + 1/4)_{j, k} = row^n_{j, k} - 
            dt/(4dx) (u^n{j + 1, k} - u^n_{j - 1, k})
        u^(n + 1/4)_{j, k} = u^n_{j, k} - 
            dt/(4dx) (p^n{j + 1, k} - p^n_{j - 1, k})
        v^(n + 1/4)_{j, k} = v^n_{j, k}
        p^(n + 1/4)_{j, k} = p^n_{j, k} - 
            (gamma * dt)/(4dx) (u^n{j + 1, k} - u^n_{j - 1, k})
    Each row is going to be represented by 1 block.

    Inputs:
        double dt_4_dN:  dt / (4dN)
        double gamma:  gamma
        int N:  N (grid size)
        double* rowOld:  device pointer to old row values
        double* uOld:  device pointer to old u values 
        double* vOld:  device pointer to old v values
        double* pOld:  device pointer to old p values

    Outputs:
        double* rowNew:  device pointer to updated row values
        double* uNew:  device pointer to updated u values 
        double* vNew:  device pointer to updated v values
        double* pNew:  device pointer to updated p values
*/
__global__ void explicitX(const double dt_4_dN, const double gamma, const int N,
    double* rowOld, double* uOld, double* vOld, double* pOld,
    double* rowNew, double* uNew, double* vNew, double* pNew)
{
    // REDACTED: kernel code is available on request
}

/*
    This is a kernel for the explicit update in y which takes the 
    following form:
        row^(n + 3/4)_{j, k} = row^(n + 1/2)_{j, k} - 
            dt/(4dx) (v^(n + 1/2)_{j, k + 1} - v^(n + 1/2)_{j, k - 1})
        u^(n + 3/4)_{j, k} = u^(n + 1/2)_{j, k}
        v^(n + 3/4)_{j, k} = v^(n + 1/2)_{j, k} - 
            dt/(4dx) (p^(n + 1/2)_{j, k + 1} - p^(n + 1/2)_{j, k - 1})
        p^(n + 3/4)_{j, k} = p^(n + 1/2)_{j, k} - 
            (gamma * dt)/(4dx) (v^(n + 1/2)_{j, k + 1} - v^(n + 1/2)_{j, k - 1})
    Each row is going to be represented by 1 block

    Inputs:
        double dt_4_dx:  dt / (4dx)
        double gamma:  gamma
        int N:  N (grid size)
        double* rowOld:  device pointer to old row values
        double* uOld:  device pointer to old u values 
        double* vOld:  device pointer to old v values
        double* pOld:  device pointer to old p values

    Outputs:
        double* rowNew:  device pointer to updated row values
        double* uNew:  device pointer to updated u values 
        double* vNew:  device pointer to updated v values
        double* pNew:  device pointer to updated p values
*/
__global__ void explicitY(const double dt_4_dN, const double gamma, const int N,
    double* rowOld, double* uOld, double* vOld, double* pOld,
    double* rowNew, double* uNew, double* vNew, double* pNew)
{
    // REDACTED: kernel code is available on request
}

/*
    This is a kernel to create the right hand side array from the four sub
    arrays.  The data passed to it is already transposed so the transposed
    right hand side can be constructed

    Inputs:
        int N:  N (number of rows)
        double* row:  device pointer to row values
        double* u:  device pointer to u values 
        double* v:  device pointer to v values
        double* p:  device pointer to p values

    Outputs:
        double* rightHandSideT:  device pointer to transposed rightHandSide
            (N rows x 4N columns)
*/
__global__ void constructRHS(int N, double* row, double* u, double* v, 
    double* p, double* rightHandSideT)
{
    // REDACTED: kernel code is available on request
}

/*
    This is a kernel to split the right hand side array into the four sub
    arrays.  The right hand side is transposed so the sub arrays will be 
    transposed.  Each block is going to handle a row in the sub arrays, 
    it will be of length N.

    Inputs:
        int N:  N (number of rows)
        double* rightHandSideT:  device pointer to transposed RHS
            (N rows x 4N columns)

    Outputs:
        double* row:  device pointer to row values
        double* u:  device pointer to u values 
        double* v:  device pointer to v values
        double* p:  device pointer to p values
*/
__global__ void splitRHS(int N, double* rightHandSideT, double* row, 
    double* u, double* v, double* p)
{
    // REDACTED: kernel code is available on request
}

int main(int argc, char* argv[])
{
    // The first process will broadcast the parameters but each process needs
    // to have memory allocated.  The problem has 2 parameters or arguments
    //      N = grid size
    //      M = time steps
    int N;
    int M;

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

    // Set the other default values
    double gamma = 1.4;
    double twoGamma = 2.0 / 1.4;
    double dt = 2.0 / M;
    double dN = 2.0 / N;
    double dt_4_dN = dt / (4 * dN);
    int fourN = 4 * N;

    // Prepare the files for saving the results
    FILE* fileid_row = NULL;
    FILE* fileid_u = NULL;
    FILE* fileid_v = NULL;
    FILE* fileid_p = NULL;

    // Open the files that will store the data
    fileid_row = fopen("EulerR.out", "wb");
    fileid_u = fopen("EulerU.out", "wb");
    fileid_v = fopen("EulerV.out", "wb");
    fileid_p = fopen("EulerP.out", "wb");

    // Initialize the timers
    struct timeval start;
    struct timeval end;

    // Start the timer
    gettimeofday(&start, NULL);

    // REDACTED: framework and implicit array initialization are available on request

    // Create the three streams, one for calling kernels, one for copying
    // implicitX, and one for calling implicitY
    cudaStream_t kernelStream;
    cudaStream_t copyXStream;
    cudaStream_t copyYStream;
    cudaStreamCreate(&kernelStream);
    cudaStreamCreate(&copyXStream);
    cudaStreamCreate(&copyYStream);

    // Start the asynchronous memory transfers
    cudaMemcpyAsync(dev_implicitX, host_implicitX, implicitSize, 
        cudaMemcpyHostToDevice, copyXStream);
    cudaMemcpyAsync(dev_implicitY, host_implicitY, implicitSize, 
        cudaMemcpyHostToDevice, copyYStream);

    // // The next step is to allocate all of the memory that will be needed on the
    // // host and device
    int N_squared = N * N;
    int size = N * N * sizeof(double);

    // Host solution arrays
    double (*host_row)[N] = (double(*)[N])malloc((N) * sizeof(*host_row));
    double (*host_u)[N] = (double(*)[N])malloc((N) * sizeof(*host_u));
    double (*host_v)[N] = (double(*)[N])malloc((N) * sizeof(*host_v));
    double (*host_p)[N] = (double(*)[N])malloc((N) * sizeof(*host_p));

    // Old arrays
    double* dev_rowOld, *dev_uOld, *dev_vOld, *dev_pOld;
    cudaMalloc(&dev_rowOld, size);
    cudaMalloc(&dev_uOld, size);
    cudaMalloc(&dev_vOld, size);
    cudaMalloc(&dev_pOld, size);

    // New arrays (intermediate results)
    double* dev_rowNew, *dev_uNew, *dev_vNew, *dev_pNew;
    cudaMalloc(&dev_rowNew, size);
    cudaMalloc(&dev_uNew, size);
    cudaMalloc(&dev_vNew, size);
    cudaMalloc(&dev_pNew, size);

    // Right hand side
    double* dev_rightHandSide;
    cudaMalloc(&dev_rightHandSide, 4 * N * N * sizeof(double));

    // Transposed implicit arrays
    double* dev_implicitXT, *dev_implicitYT;
    cudaMalloc(&dev_implicitXT, implicitSize);
    cudaMalloc(&dev_implicitYT, implicitSize);

    cudaMemcpyToSymbol(dev_N, &N, sizeof(int));
    cudaMemcpyToSymbol(dev_gamma, &gamma, sizeof(double));
    cudaMemcpyToSymbol(dev_twoGamma, &twoGamma, sizeof(double));
    cudaMemcpyToSymbol(dev_dN, &dN, sizeof(double));
    cudaMemcpyToSymbol(dev_dt_4_dN, &dt_4_dN, sizeof(double));


    // // Call the kernel to set the initial conditions
    initialize<<<N, N>>>(dev_rowOld, dev_uOld, dev_vOld, dev_pOld);

    // Create streams for the ADI solver
    cudaStream_t stream1;
    cudaStream_t stream2;
    cudaStream_t stream3;
    cudaStream_t stream4;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);
    cudaStreamCreate(&stream4);

    // Synchronize and free the streams
    cudaStreamSynchronize(kernelStream);
    cudaStreamSynchronize(copyXStream);
    cudaStreamSynchronize(copyYStream);

    cudaStreamDestroy(kernelStream);
    cudaStreamDestroy(copyXStream);
    cudaStreamDestroy(copyYStream);

    // Transpose implicit arrays (4N x 4N)
    transpose<<<gridImpl, blockImpl, 0, stream1>>>(dev_implicitX, fourN, dev_implicitXT);
    transpose<<<gridImpl, blockImpl, 0, stream2>>>(dev_implicitY, fourN, dev_implicitYT);

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);
    cudaStreamSynchronize(stream4);

    // cuSolverDn is going to be used to replicate the LU decomp and solutions
    // utilizing the GPU.
    
    // REDACTED: cuSolver and getrf code are available on request

    // Copy the initial conditions from the device to the host and save
    // the data
    cudaMemcpy(host_row, dev_rowOld, N_squared * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_u, dev_uOld, N_squared * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_v, dev_vOld, N_squared * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_p, dev_pOld, N_squared * sizeof(double), cudaMemcpyDeviceToHost);

    // Save the initial data
    fwrite(host_row, sizeof(double), (N * N), fileid_row);
    fwrite(host_u, sizeof(double), (N * N), fileid_u);
    fwrite(host_v, sizeof(double), (N * N), fileid_v);
    fwrite(host_p, sizeof(double), (N * N), fileid_p);

    // dim3 gridSol(N / tile_dim, N / tile_dim);

    // // // ADI Algorithm
    // // Information about the arrays being used (number of rows x number of columns)
    //      rowOld/uOld/vOld/pOld = old solution (N x N)
    //      rowNew/uNew/vNew/pNew = new solution (N x N)
    //      rightHandSide = right hand side of implicit steps (4N x localRows)
    //      rightHandSideT = transpose of RHS (localRows x 4N)
    // For this problem there are periodic boundary conditions so the edge cases
    // are going to look like this where X is one of {row, u, v, p}
    //      X^n_{-1, k} = X^n_{N - 1, k}
    //      X^n_{N, k} = X^n_{0, k}
    //      X^n_{j, -1} = X^n_{j, N - 1}
    //      X^n_{j, N} = X^n_{j, 0}

    // REDACTED: ADI algorithm is available on request

    // End the timer and then determine the max time
    gettimeofday(&end, NULL);

    double elapsedTime = (end.tv_sec + end.tv_usec / 1.0e6) - 
        (start.tv_sec + start.tv_usec / 1.0e6);

    printf("Elapsed time: %.6f seconds\n", elapsedTime);

    // Free cuSolverDn data
    cudaFree(dev_work);
    cudaFree(dev_ipivX);
    cudaFree(dev_ipivY);
    cudaFree(dev_info);
    cusolverDnDestroy(cusolverH);

    // Close the files, free the memory and save the results
    fclose(fileid_row);
    fclose(fileid_u);
    fclose(fileid_v);
    fclose(fileid_p);

    // Print parameters and results
    double timePerIteration = elapsedTime / M;
    printf("Grid Length: %d\n", N);
    printf("Time Steps: %d\n", M);
    printf("Function took %f seconds to execute\n", elapsedTime);
    printf("Average Time per Iteration: %lf\n", timePerIteration);

    // Free the data used
    cudaFree(host_implicitX);
    cudaFree(host_implicitY);
    free(host_implicitXT);
    free(host_implicitYT);
    cudaFree(dev_implicitX);
    cudaFree(dev_implicitY);
    free(host_row);
    free(host_u);
    free(host_v);
    free(host_p);
    cudaFree(dev_rowOld);
    cudaFree(dev_uOld);
    cudaFree(dev_vOld);
    cudaFree(dev_pOld);
    cudaFree(dev_rowNew);
    cudaFree(dev_uNew);
    cudaFree(dev_vNew);
    cudaFree(dev_pNew);
    cudaFree(dev_rightHandSide);
    cudaFree(dev_implicitXT);
    cudaFree(dev_implicitYT);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    cudaStreamDestroy(stream4);

    return 0;
}