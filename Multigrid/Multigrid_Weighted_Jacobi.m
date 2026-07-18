% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, total_cost, error] = Multigrid_Weighted_Jacobi(u, N, h, ...
    omega, delta, b, n, iter_cost, num_iterations)
    % This function uses weighted Jacobi to solve a system of linear
    % equations which is generalized for systems of any size.  This is a
    % one dimensional problem so there is no need to worry about the
    % mapping of u.  In addition, the sparsity and discretization of A
    % allows for more efficient computing in the solver and matrix
    % multiplications the latter of which being handled in a separate
    % function.  For the weighted Jacobi implementation there will be a
    % vector of the old and new u values
    % Inputs:
    %   u = the initial data for the system
    %   N = the number of data points in u, does include the boundary
    %   points
    %   h = step size in the grid
    %   omega = weight term
    %   delta = parameter from original boundary value problem
    %   b = solution to the boundary value problem
    %   n = parameter from original problem
    %   iter_cost = computational cost per Jacobi pass
    %   num_iterations = number of Jacobi steps to take
    % Outputs
    %   u_final = the results of the relaxation method for this grid level
    %   total_cost = computational cost for this Jacobi solver

    % This is the form of weighted Jacobi:
    %   x_(n+1) = (1-w)x_n + w(Gx_n + D^-1b)
    % This method converges when the RMS norm of U_n - U_exact is less than
    % the tolerance.  Since the U_exact vector depends on the size of the
    % grid it will be calculated and used at each step of the method. The
    % exact u_j = C_tilda(nh) * sin(n * pi x_j) where x_j = jh and c_tilda 
    % are determined from the problem.  U_exact is going to be the same
    % size as u and will be populated below.
    % to be the same size as u and I am going to populate the values here:
    tol = 10^-5;
    error = 1;
    u_exact = u;
    u_exact(1) = 0;
    u_exact(end) = 0;
    % These are used for the c_tilda calculations
    first_term = n^2 * pi^2;
    second_term = (sin(1/2 * n * pi * h))^2;
    third_term = (1/4) * n^2 * pi^2 * h^2;
    c_tilda = 1/(first_term * (second_term / third_term) + delta);
    for i = 2:(N - 1)
        u_exact(i) = c_tilda * sin(n * pi * (i - 1) * h);
    end

    % We need to keep 2 sets of data, the old points and the new points
    u_old = u;
    b;
    u_new = u;

    % For weighted Jacobi, and relaxation methods in general, A is the same
    % for each iteration.  In this case the structure is defined so we can
    % precalculate some of the values.  Here are the steps I am going to
    % use for weighted Jacobi:
    % These are done outside of the main loop, the values do not change
    %   1) G = I - D^-1 * A
    %   2) G(w) = (1 - w) * I + w * G
    %   3) D_b = D^-1 * b
    % These are done inside of the while loop for each value of u
    %   4) u_(n+1) = (1 - w) * u_n + w(G * u_n + D_b)
    % This process is going to be run until the relaxation method
    % converges, or when the maximum value of u is less than the tolerance
    % value.

    % From the discretization of A, this is what each row of A looks like:
    %   [ ... 0 0 -1 (2 + h^2*delta^2) -1 0 0 ... ]
    % As a result each row of G is going to look like:;
    %   [ ... 0 0 1 / (2 + h^2*delta^2) 0 1 / (2 + h^2*delta^2) 0 0 ... ]
    % From here we can determine what each row of G(w) will look like,
    %   [ ... 0 0 w / (2 + h^2*delta^2) (1 - w) w / (2 + h^2*delta^2) 0 0
    %   ... ]

    % The D_b vector is easy to calculate as each value is 
    % b / (2 + h^2*delta^2)
    D_b = b;
    for i = 1:N
        D_b(i) = b(i) / (2 + h^2 * delta^2);
    end
    D_b;

    % REDACTED: the weighted Jacobi code is available on request

    error;
    % Once the iteration is complete we determined the u values that lead
    % to convergence so we can return them
    u_final = u_new;
    u_exact;
end