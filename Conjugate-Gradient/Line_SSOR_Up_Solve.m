% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_gs] = Line_SSOR_Up_Solve(D, p, omega, u_row, prior_u_row, b, num_grid_points)
    % This function is going to be similar to a Thomas algorithm except for
    % solving a system with 2 sub diagonals.  This problem models the up
    % solve for the line SSOR preconditioner.  The system being solved is
    % of the form 
    %   Au_gs + u_New_(i+1) = b
    % In this case the values of u_New are known and the coefficient is 1,
    % so the problem being solved is, 
    %   Au_gs = b - u_New_(i+1)
    % where A is a diagonal and a super diagonal so it requires 1 sweep to
    % solve the system.  One assumption being made about this problem is
    % that the boundary conditions can be ignored since they are all 0, a
    % row of A looks like,
    %       [0 ... 0 0 D wp ... 0]
    % Inputs:
    %   D = diagonal value
    %   p = parameter from original PDE
    %   omega = SOR weight term
    %   u_row = old values of u for the current row
    %   prior_u_row = new u values from the prior row
    %   b = parameter from original equation
    %   num_grid_points = equations per row
    % Outputs:
    %   u_gs = solution to the system

    % The first step is to initialize the return value and making sure that
    % the end conditions are satisfied.  The rest of the function is going
    % to use the structure of A to solve for the solution to the system.

    u_gs = u_row;
    u_gs(1) = 0;
    u_gs(end) = 0;

    % REDACTED: the rest of the line SSOR up solve code is available on
    % request
end