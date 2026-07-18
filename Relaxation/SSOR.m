% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = SSOR(u, epsilon, N, omega)
    % This function is going to use SSOR to solve a system of linear 
    % equations.  The sparsity of A will avoid a mapping of the 2D A to a 
    % 1D vector and an optimized matrix multiplication using the structure 
    % of A
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

    % SSOR converges when the largest u value is below the tolerance, with 
    % the starting value being set to the initial max.  The total number of 
    % equations used in the solver is N^2 and the system dimension is N + 2 
    % when accounting for the ghost points
    tol = 10^-7;
    max_val = max(max(u));
    num_eqn = N^2;
    sys_dim = N + 2;

    % Three sets of data are stored, the initial values, the intermediate
    % step, and the final step with all three being the same size, (N + 2) 
    % x (N + 2) = N_mod x N_mod
    u_initial = u;
    u_intermediate = u;
    u_full_step = u;

    % SSOR is built off of SOR, using a SOR sweep down and then a reverse
    % sweep of SOR up.  The first sweep generates a temporary set of u
    % values which will be used in the second sweep to generate the final u
    % values for the step.  This program consists of 2 nested for loops for
    % the down and up sweep with each having N^2 operation.  On the up
    % sweep the indices will be reversed but the structure will stay the
    % same and both sweeps will take the form:
    %   u_i^(n+1) = omega * (b_i - sum(not including i) of a_ij*x_j) / a_ii
    %   + u_i^n
    % The while loop will run until the method converges or the maximum
    % number of iterations are completed
    iterations = 0;

    % REDACTED: SSOR code is available on request

    % After the loop the final u values are found and returned
    u_final = u_new; 
end