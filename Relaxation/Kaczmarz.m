% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, num_iterations, max_val] = Kaczmarz(u, epsilon, N, ...
    iterations)
    % This function uses the Kaczmarz method to solve a system of linear
    % equations 
    % Inputs:
    %   u = initial data for the system
    %   epsilon = constant from original PDE
    %   N = number of equations
    %   iterations = maximum iterations to run solver for
    % Outputs:
    %   u_final = the results of the relaxation method
    %   num_iteration = number of iterations solver ran for
    %   max_val = ending maximum value of u 

    % In the Kaczmarz method, each step of the iteration involves updating
    % all of the u values resulting in N^2 updates to u.  The code is going
    % to be contained in a single while loop that checks for convergence
    % and number of steps in the iterations.  This is the general structure
    % of the loop:
    %   for k = 1:iterations????
    %       for i = 1:N
    %           r_i = b_i - dot_product(a_i, x_k
    %           delta_i= r_i / dot_product(a_i, a_i)
    %           x^k = x^k + delta_i * a_i
    %       end
    %       x^(k + 1) = x^k
    %   end
    % Using the structure of A these are the possible dot product values:
    % [-4, 1, 1] = 18
    % [-4, 1, 1, 1] = 19
    % [-4, 1, 1, 1, 1] = 20

    % Flatten u into a column vector ignoring the boundary points
    num_equations = N^2;
    flat_u = zeros(num_equations, 1);
    cur_index = 1;
    for i = 2:(N+1)
        for j = 2:(N + 1)
            flat_u(cur_index) = u(i, j);
        end
    end

    % These are values used in the loop
    max_val = max(flat_u);
    max_val = max(max(u));
    tol = 10^-7;
    num_iterations = 0;
    residual = 0;
    delta_i = 0;
    diag_coef = -2 - 2 * epsilon;

    % REDACTED: the Kaczmarz loop is available on request

    u_final = u;
end
