% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [y] = SSOR_Preconditioner(omega, b, N)
    % This function is going to implement SSOR to solve for the
    % preconditioning step.  There is a single forward sweep and a
    % backwards sweep that follows it.  The first sweep is solving 
    %   (D - omega * L)x = omega * b
    % This is going to be solved using SOR, and the backwards sweep is
    % going to be of the form
    %   (D - omega * L^T)y = omega * b + omega * L * x + (1 - omega)Dx
    % in this case y is what we are solving for and it the value of M^-1b.
    % Also for this case A = D + L + L^T.  
    % Inputs:
    %   omega = SOR weight term
    %   b = solution to the system, includes boundary points
    %   N = number of equations, includes the boundary points
    % Outputs:
    %   y = the value of M^-1b

    % From the discretization of 4 we have the following for the rows of A
    %       [0 ... 0 -1 2 -1 0 ... 0]
    % This means that D is 2I, L is the sub diagonal, and L^T is the super
    % diagonal.  

    % The initial guess for x is going to be all zeros and this is the
    % first sweep to solve (D - omega * L)x = omega * b.
    x_initial = zeros(N, 1);

    % REDACTED: SSOR preconditioner code is available on request

end