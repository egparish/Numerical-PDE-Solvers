% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, error, computational_cost] = V_Cycle(b, start_u, alpha, ...
    num_grids, delta, n, v_down, v_up, epsilon)
    % This function is going to run the V_cycle given an initial condition
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

    % There will be three vectors that will store all of the needed
    % information, the Q vector stores the x values, the F vector stores the
    % b values, and a start_indices vector that stores the index where each
    % of the grid values start in the Q and F vectors.  There will be two
    % additional vectors, size which store the size of each grid, and
    % end_indices, which stores the ending index of each subset of F and Q.
    % All of the vectors will start with the finest grid having the first
    % values and working down to the coarser grids.  The sizes vector is
    % going to store the (N + 2) values which gives N interior points and 2
    % boundary points.  The start_indices vector is going to store the true
    % starting indices, since MATLAB indexing starts at 0 this will be the
    % first value and then the grid sizes will be added to it.  The
    % end_indices vector is going to store the MATLAB index of the last
    % value which is one subtracted from the size added to the starting
    % index.

    % The first step is to initialize and fill the vectors that will be used
    % which depends on alpha and the number of grids that will be used.
    % For this implementation, alpha is going to be the number of internal
    % points + 1 in the coarsest grid.  The rest of the constants used in
    % the solver are also defined
    % alpha = N + 1 of smallest grid size, N = number of interior points
    % % % alpha = 4;
    % % % num_grids = 6;
    % % % delta = 1;
    % % % n = 4;
    % % epsilon = 10^-5;
    omega = 0.8;
    max_grid_size = alpha * 2^(num_grids - 1) + 1;
    h_start = 1 / (max_grid_size - 1);
    total_cost = 0;
    % % v_down = 100;
    % % v_up = 100;

    % Fill the vectors remembering the finest grid is first, with
    % num_interior_points storing the number of interior points for each
    % grid
    sizes = zeros(1, num_grids);
    start_indicies = zeros(1, num_grids);
    end_indicies = zeros(1, num_grids);
    h_values = zeros(1, num_grids);
    num_interior_points = zeros(1, num_grids);
    cost_per_sweep = zeros(1, num_grids);

    % The total length is the number of all internal and boundary points
    % across all of the grids.  Each grid has 2 boundary points and the set
    % of internal points being stored.  While the boundary points are not
    % used for this problem it is easier to store them.  For each grid, the
    % following will be used to calculate the total number of points at
    % that level, N + 1 = alpha * 2^(grid_num - 1) + 1
    total_length = 0;
    
    % Populate the values using a loop
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

    % The first part of multigrid involves passing down through the
    % V-cycle, update the current x value and the next b values.  This
    % starts with the finest grid resulting in (num_grids - 1) iteration,
    % with each iteration accessing the current b and x values along with
    % the future b values.  The last x value will be updated after the down
    % loop.
    grid_num = 1;

    % REDACTED: down V-cycle is available on request
    
    % After reaching the bottom of the cycle, gaussian elimination is used
    % on the coarsest grid to solve for the last x values
    F_final = F(start_indicies(end):end_indicies(end));
    X_final = Multigrid_Gaussian_Elimination(F_final, sizes(end), h_values(end), delta);
    Q(start_indicies(end):end_indicies(end)) = X_final;

    % After the Gaussian elimination the V-cycle can be passed through
    % moving upwards.  These iterations start with the coarsest grid
    % working on the current Q and F values along with the next Q values to
    % update them.  After the cycle is done it is zeros out as it is no
    % longer needed.  The final u values are found arriving to the finest
    % grid with nothing being done the second time the finest grid is
    % visited

    % REDACTED: up V-cycle is available on request

    computational_cost = total_cost;
    u_final = current_x;
    error = error_up;
end