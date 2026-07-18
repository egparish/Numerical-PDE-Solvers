% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, iterations, max_val] = Line_SSOR(u, epsilon, N, omega)
    % This function uses the line SSOR method to solve a system of linear
    % equations.  Since A is derived from the discretization of the problem
    % it is not stored and a separate function is used for matrix
    % multiplications involving A
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

    % Line SSOR is similar to line SOR but it includes an additional sweep
    % up after the initial sweep down.  Line SOR involves solving all of the
    % equations in a line simultaneously using the Thomas algorithm, as a
    % result the Thomas algorithm will be used in each step to solve for
    % the entire row of u.  This method converges when the largest u values
    % is below the tolerance with the starting value being set to the
    % initial max.  The total number of equations being used in the solver
    % is N^2, when including the ghost points the total system dimension is 
    % N + 2
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

    % As mentioned above, line SSOR is built off of line SOR, involving a
    % line SOR sweep down and then a reverse line sweep of SOR up.  The
    % first sweep generates a set of temporary u value which will be used
    % in the second sweep to generate the final u values for the step.
    % This program is going to consist of 2 nested for loops, with each
    % calling the Thomas algorithm on the line of u it is working with.
    % the main while loop will run until the method converges or the
    % maximum number of iterations are completed.  The first sweep down is
    % the same as line SOR with the sweep up involving the indices be
    % reversed.  Both of the sweeps are going to use the form:
    %   GS_vals = Thomas_Algorithm(current row)
    %   u_new(i, row) = omega * GS_vals(i) + (1 - omega) * u_old(i, row)
    iterations = 0;

    % REDACTED: line SSOR code is available on request

    % After the loop the final u values are found and returned
    u_final = u_initial;
end