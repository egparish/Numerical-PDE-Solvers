% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, total_cost, error] = Multigrid_Gauss_Seidel(u, N, h, ...
    delta, b, n, iter_cost, num_iterations)
    % This function is going to use Gauss-Seidel to solve the system but
    % one major difference is that it will work with any sized system and
    % go off of one the passed parameters.  This is a one dimensional
    % problem so there is no concern about the mapping of u and the
    % discretization of the original problem makes the system easier to
    % work with.
    % Inputs:
    %   u = the initial data for the system
    %   N = the number of data points in u, does include the boundary
    %   points
    %   h = step size in the grid
    %   delta = parameter from original boundary value problem
    %   b = solution to the boundary value problem
    %   n = parameter from original problem
    %   iter_cost = computational cost per Jacobi pass
    %   num_iterations = number of Jacobi steps to take
    % Outputs
    %   u_final = the results of the relaxation method for this grid level
    %   total_cost = computational cost for this Gauss-Seidel solver

    % This function will use one set of u values and the new values will be
    % used as soon as they are determined.  Within the Gauss-Seidel loop
    % the interior points will be iterated through which represent the
    % equations being solved.  Each new u value will be found as follows:
    %   u_i^(n+1) = (b_i - sum(not including i) of a_ij*x_j) / a_ii
    % The loop will run until the maximum number of iterations or until an
    % error tolerance is met.  The error is when the RMS norm of U_n -
    % U_exact is less than the tolerance.  Since the U_exact vectors
    % depends on the size of the grid it will be calculated at each step of
    % the grid.  The exact u_j = C_tilda(nh) * sin(n * pi x_j) where
    % x_j = jh and C_tilda is determined from the problem with the values
    % being found below.
    u_vals = u;

    % REDACTED: the Gauss-Seidel code is available on request
   
    % Once the iteration is complete we determined the u values that lead
    % to convergence so we can return them
    u_final = u_vals;
    u_exact;
end