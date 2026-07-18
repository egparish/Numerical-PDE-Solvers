% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [x_output] = Multigrid_Gaussian_Elimination(b, N , h, delta)
    % When we get to the bottom of the V-cycle we are going to use gaussian
    % elimination to solve Ax = b on the coarsest grid.  The structure of A
    % is known and the b values come from the final grid transfer of the
    % residual.  The Thomas algorithm can be used for Gaussian elimination
    % since A is tridiagonal.  The will give the final x values which can
    % be used to move up the V-cycle
    % Inputs:
    %   b = residual from the coarsest grid
    %   N = size of the coarsest grid, includes the boundary points
    %   h = step size in the grid
    %   delta = parameter from original boundary value problem
    % Outputs:
    %   x_output = the solution to Ax = b for the coarsest grid

    % The Thomas algorithm involves a forward and backwards sweep with a
    % scratch vector being needed.  From the discretization of A, each row
    % looks like:
    %   [ ... 0 0 -1 (2 + h^2*delta^2) -1 0 0 ... ]
    x_output = zeros(N, 1);
    scratch = zeros(N, 1);
    a = -1; 
    c = -1; 
    d = 2 + h^2 * delta^2;

    % REDACTED: the Thomas algorithm code is available on request

end