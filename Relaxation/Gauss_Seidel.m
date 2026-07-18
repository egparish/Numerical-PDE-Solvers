% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = Gauss_Seidel(u, epsilon, N)
    % This function is going to use Gauss-Seidel to solve a system of
    % linear equations.  Since A is a sparce matrix, a separate function is
    % going to be used to optimize the matrix multiplication.  In addition,
    % A does not need to be stored since it is dependent on the
    % discretization of the problem.  The u matrix is going to be
    % (N + 2) x (N + 2) = N_mod x N_mod
    % Inputs:
    %   u = the initial data for the system
    %   epsilon = from the original PDE
    %   N = the number of data points in u, does not include the imaginary
    %   points
    % Outputs
    %   u_final = the results of the relaxation method

    % This method converges when the largest u value is below the tolerance
    % which starts as the initial max value.  The total number of equations
    % used in the solver is N^2.  The system_dimension variable is the size
    % of the matrix including the fake points, sys_dim = N + 2
    tol = 10^-7;
    max_val = max(max(u));
    num_eqn = N^2;
    sys_dim = N + 2;

    % In Gauss-Seidel, the new values are used as soon as they are found
    % which requires only one matrix for the solution.  The solver will be
    % run using a while loop until the tolerance is met or the maximum
    % number of iterations have been completed.  Within the solver loop
    % there is a double loop that mimics iterating through all of the
    % equations which is done due to how u is stored.  Each new value of u
    % will be found as follows:
    %   u_i^(n+1) = (b_i - sum(not including i) of a_ij*x_j) / a_ii
    u_vals = u;
    iterations = 0;

    % REDACTED: Gauss-Seidel loop is available on request

    % After the loop the final u values were found
    u_final = u_vals;
end