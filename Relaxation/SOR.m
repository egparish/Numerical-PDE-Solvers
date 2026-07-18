% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = SOR(u, epsilon, N, omega)
    % This function is going to use SOR to solve a system of linear 
    % equations.  The sparsity of A will avoid a mapping of the 2D A to
    % a 1D vector and an optimized matrix multiplication using the
    % structure of A
    % Inputs:
    %   u = the initial data for the system
    %   epsilon = from the original PDE
    %   N = the number of data points in u, does not include the imaginary
    %   points
    %   omega = SOR weight
    % Outputs
    %   u_final = the results of the relaxation method
    %   num_iteration = number of iterations solver ran for
    %   max_val = ending maximum value of u 

    % SOR converges when the largest u value is below the tolerance, with 
    % the starting value being set to the initial max.  The total number of 
    % equations used in the solver is N^2 and the system dimension is N + 2 
    % when accounting for the ghost points
    tol = 10^-7;
    max_val = max(max(u));
    num_eqn = N^2;
    sys_dim = N + 2;

    % Since new values are used as they are generated only one set of data
    % is needed with the size of u being (N + 2) x (N + 2) = N_mod x N_mod
    u_vals = u;

    % This solver is similar to the Gauss-Seidel code with a single double
    % loop, with N^2 steps, being used to update all of the values of u.
    % The main difference with SOR is that there is a weight applied to the
    % Gauss-Seidel sum and then the old u value is added to it.  Each new u 
    % value is going to be found as follows:
    %   u_i^(n+1) = omega * (b_i - sum(not including i) of a_ij*x_j) / a_ii
    %   + u_i^n
    % The while loop will run until the method converges or the maximum
    % number of iterations are completed.
    iterations = 0;

    % REDACTED: SOR code is available on request

    % After the loop the final u values are found and returned
    u_final = u_new;

end