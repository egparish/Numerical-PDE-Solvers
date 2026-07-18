% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [b_new] = Fine_to_Coarse_Grid_Trasnsfer(r_old, N_old, N_new)
    % This function is going to handle the fine to coarse grid transfer
    % when going down the V-Cycle.  This is applied to the calculated
    % residual to find the b values for the next step.  When moving down
    % the even points are used as the base with the some weight from the
    % odd points to not lose that information.
    % Inputs:
    %   r_old = the residuals from the finer grid
    %   N_old = the size of the r_old vector, includes the boundary points
    %   N_new = the size of the new b vector, includes the boundary points
    % Outputs
    %   b_new = the result of the fine to coarse grid transfer

    % This function involves initializing the b_new vector and populating
    % it with values from the finer grid.  Since there are Dirichlet
    % boundary conditions those can be left alone to have a value of 0.
    % The next step involves iterating through the r_old vector and
    % populating the b_new vector.  Full weighting will be used to not lose
    % the information from the odd points which involves half of the weight
    % from the even point and a quarter of the weight from the adjacent
    % points.
    b_new = zeros(1, N_new);
    b_new_index = 1;

    % REDACTED: find to coarse transfer code is available on request
    
end