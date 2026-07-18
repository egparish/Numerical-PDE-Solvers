% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [x_new] = Coarse_to_Fine_Grid_Trasnsfer(x_old, N_old, N_new)
    % This function is going to handle the coarse to fine grid transfer
    % when going up the V-Cycle.  This is applied to old x values to find
    % the x values for the finer grid.  When moving up we need to create
    % new values through extrapolation.
    % Inputs:
    %   x_old = the x values from the coarse grid
    %   N_old = the size of the x_old vector, includes the boundary points
    %   N_new = the size of the new x vector, includes the boundary points
    % Outputs
    %   x_new = the result of the coarse to fine grid transfer

    % This function will initialize the x_new vector, and populate the new
    % vector using the old x values.  All of the boundary points and even
    % indexed points will directly transfer to the new vector.  For the odd
    % indexed points it will be the average of the surrounding points.  The
    % Dirichlet conditions will be handled outside of the loop and will be
    % set to 0.
    x_new = zeros(1, N_new);
    next_x_to_move = 1;

    % REDACTED: coarse to fine grid transfer code is available on request

end