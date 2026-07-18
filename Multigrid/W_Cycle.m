% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, error, computational_cost] = W_Cycle(b, start_u, alpha, ...
    num_grids, delta, n, v_down, v_up, epsilon)
    % This function is going to run a W-cycle given an initial condition.
    % A W-cycle is when a V-cycle is run each time an intermediate grid is
    % visited for the first time when moving up.  This is tracked by
    % keeping flag variables for each layer
    % Inputs:
    %   b = the solution to the boundary value problem, this is used for
    %   convergence checks
    %   start_u = the starting guess for u
    %   alpha = N + 1 of smallest grid size, N = number of interior points
    %   num_grids = the number of grids the V-Cycle is run on
    %   delta = parameter from the original problem
    %   n = parameter from the original problem
    %   v_down = number of relaxation iterations to run when moving down
    %   the V-Cycle
    %   v_up = number of relaxation iterations to run when moving up the
    %   V-Cycle
    %   epsilon = tolerance check value
    % Outputs:
    %   u_final = the final u values
    %   error = the error from one V-Cycle
    %   computational_cost = computational cost from one V-cycle

    % This implementation is going to be modified from the V-cycle
    % function, with the changes being a flag vector and a variable that
    % tracks the current grid where 1 is the finest grid and num_grids is
    % the coarsest grid.  The W-cycle is finished when the coarsest grid is
    % visited for the second time which will be tracked with a counter and
    % the while loop.  When the coarsest grid is visited, the only option
    % is to move up and the flag value is irrelevant.  The starting
    % initialization is used from the V-cycle function
    max_grid_size = alpha * 2^(num_grids - 1) + 1;
    h_start = 1 / (max_grid_size - 1);
    total_cost = 0;
    omega = 0.8;
    sizes = zeros(1, num_grids);
    start_indicies = zeros(1, num_grids);
    end_indicies = zeros(1, num_grids);
    h_values = zeros(1, num_grids);
    num_interior_points = zeros(1, num_grids);
    cost_per_sweep = zeros(1, num_grids);
    total_length = 0;

    % This populates all of the vectors we are going to use
    for i = 1:num_grids
        % The sizes vector is easy to determine
        sizes(i) = alpha * 2^(num_grids - i) + 1;
        total_length = total_length + sizes(i);
        num_interior_points(i) = sizes(i) - 2;
        h_values(i) = 1 / (sizes(i) - 1);
        cost_per_sweep(i) = 1 * 2^(-i + 1);
    
        % The indices vector has edge cases
        if i == 1
            start_indicies(i) = 1;
            end_indicies(i) = start_indicies(i) + sizes(i) - 1;
        else
            start_indicies(i) = end_indicies(i - 1) + 1;
            end_indicies(i) = start_indicies(i) + sizes(i) - 1;
        end
    end

    % Afterwards the F and Q vectors can be initialized to be all zeros
    F = zeros(1, total_length);
    Q = zeros(1, total_length);

    % Using the data structures for F and Q described above they can be
    % populated with the b values for the interior points of F and the
    % initial guess for Q for the finest grid
    F(start_indicies(1):end_indicies(1)) = b;
    Q(start_indicies(1):end_indicies(1)) = start_u;

    % This is the initialization specific to the W-Cycle
    finest_grid_visited = 0;
    current_grid = 1;
    flags = ones(1, num_grids);
    num_grids;

    % REDACTED: W-cycle code is available on request

    % After the while loop finished the return values need to be updates
    computational_cost = total_cost;
end