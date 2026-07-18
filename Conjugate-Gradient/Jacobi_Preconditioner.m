function [y] = Jacobi_Preconditioner(b, N)
    % This function is going to implement Jacobi, or diagonal scaling, as
    % the preconditioner.  It involves using the diagonals of A to scale
    % the residues
    % Inputs:
    %   b = solution to the system, includes boundary points
    %   N = number of equations, includes the boundary points
    % Outputs:
    %   y = the value of M^-1b

    % From the discretization of 4 we have the following for the rows of A
    %       [0 ... 0 -1 2 -1 0 ... 0]
    % This means that each value in b is going to be divided by 2

    y = b;

    for i = 1:N
        y(i) = y(i) / 2;
    end
end