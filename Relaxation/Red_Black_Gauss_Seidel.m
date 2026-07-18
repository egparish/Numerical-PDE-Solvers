% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = Red_Black_Gauss_Seidel(u, epsilon, N)
    % This function is going to use red-black Gauss-Seidel to solve a
    % system of linear equations.  The sparsity of A will avoid a mapping 
    % of the 2D A to a 1D vector and an optimized matrix multiplication 
    % using the structure of A.  For the red-black indexing, all of the
    % would be odd indices without the ghost points are going to be red and
    % the would be even indices are going to be black
    % Inputs:
    %   u = the initial data for the system
    %   epsilon = from the original PDE
    %   N = the number of data points in u, does not include the imaginary
    %   points
    % Outputs
    %   u_final = the results of the relaxation method
    %   num_iteration = number of iterations solver ran for
    %   max_val = ending maximum value of u 

    % Red-Black Gauss-Seidel converges when the largest u value is below the
    % tolerance, with the starting value being set to the initial max.  The
    % total number of equations used in the solver is N^2 and the system
    % dimension is N + 2 when accounting for the ghost points
    tol = 10^-7;
    max_val = max(max(u));
    num_eqn = N^2;
    sys_dim = N + 2;

    % Since new values are used as they are generated only one set of data
    % is needed with the size of u being (N + 2) x (N + 2) = N_mod x N_mod
    u_vals = u;

    % This solver is similar to the Gauss-Seidel code with the main
    % difference is that there are two sets of nested double loops.  The
    % total number of computations is sill N^2 but the solver is split into
    % red updates and then black updates with the new values being used as
    % soon as they are found.  Each new  u value is going to be found as 
    % follows:
    %   u_i^(n+1) = (b_i - sum(not including i) of a_ij*x_j) / a_ii
    % The while loop will run until the method converges or the maximum 
    % number of iterations are completed. 
    iterations = 0;

    % REDACTED: red-black Gauss-Seidel code is available on request

    % After the loop the final u values are found and returned
    u_final = u_new; 
end