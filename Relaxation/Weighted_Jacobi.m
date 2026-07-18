% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = Weighted_Jacobi(u, epsilon, N, omega)
    % This function is going to use weighted Jacobi to solve a system of
    % linear equations.  The sparsity of A will avoid a mapping of the 2D A
    % to a 1D vector and an optimized matrix multiplication using the
    % structure of A
    % Inputs:
    %   u = the initial data for the system
    %   epsilon = from the original PDE
    %   N = the number of data points in u, does not include the imaginary
    %   points
    %   omega = weight term
    % Outputs
    %   u_final = the results of the relaxation method
    %   num_iteration = number of iterations solver ran for
    %   max_val = ending maximum value of u 

    % Weighted Jacobi converges when the largest u value is below the
    % tolerance, with the starting value being set to the initial max.  The
    % total number of equations used in the solver is N^2 and the system
    % dimension is N + 2 when accounting for the ghost points
    tol = 10^-7;
    max_val = max(max(u));
    num_eqn = N^2;
    sys_dim = N + 2;

    % Two sets of data are stored, the old points and the new points.  Both
    % are the same size of u and are (N + 2) x (N + 2) = N_mod x N_mod
    u_old = u;
    u_new = zeros(sys_dim, sys_dim);

    % Weighted Jacobi is similar to point Jacobi with the updates taking
    % the form:
    %   x_(n+1) = (1-w)x_n + w(Gx_n + D^-1b)
    % The while loop will run until the method converges or the maximum
    % number of iterations are completed
    iterations = 0;

    % REDACTED: weighted Jacobi code is available on request

    % After the loop the final u values are found and returned
    u_final = u_new;   
end