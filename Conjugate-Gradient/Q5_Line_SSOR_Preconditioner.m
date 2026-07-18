% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [y] = Q5_Line_SSOR_Preconditioner(omega, b, N, p, num_grid_points)
    % This function is going to implement line SSOR to solve for the
    % preconditioning step for question 5.  There is a single forward 
    % sweep and a backwards sweep that follows it.  The first sweep is 
    % solving 
    %   (D - omega * L)x = omega * b
    % This is going to be solved using line SOR, and the backwards sweep is
    % going to be of the form
    %   (D - omega * L^T)y = omega * b + omega * L * x + (1 - omega)Dx
    % in this case y is what we are solving for and it the value of M^-1b.
    % Also for this case A = D + L + L^T.  
    % Inputs:
    %   omega = SOR weight term
    %   b = solution to the system, includes boundary points
    %   N = number of equations, includes the boundary points
    %   p = parameter from original equation
    %   num_grid_points = equations per row
    % Outputs:
    %   y = the value of M^-1b

    % From the discretization of 3 we have the following for the rows of A
    %            [0 ... -1 ... 0 -p 2(p + 1) -p 0 ... -1 ... 0]
    % This means that D is 2(p + 1), L is composed of the sub diagonal of
    % -p and another sub diagonal.  L^T is composed of the super diagonal
    % of -p and another super diagonal

    % The initial guess for x is going to be all zeros and this is the
    % first sweep to solve (D - omega * L)x = omega * b.  
    x_initial = zeros(N, 1);
    x_intermediate = x_initial;
    y = x_initial;

    % REDACTED: line SSOR preconditioner code is available on request

end