% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [A_x_output] = Matrix_Multiplication(x_input, delta, h, N)
    % This function is going to handle the sparce matrix multiplication
    % associated with the boundary value problem.  The structure of A is
    % known and the values are constant so A does not need to be stored.
    % The interior points can all be treated the same but the boundary
    % conditions are different so they will involve 2 terms instead of 3.
    % Inputs:
    %   x_input = the x vector we want to multiply with A
    %   delta = the value of the delta parameter from the original question
    %   h = spacing of the grid
    %   N = size of x
    % Outputs:
    %   A_x_output = the result of A * x, is the same size as x

    % The structure of the original boundary value problem is:
    %   -u_xx + delta^2*u = sin(n*pi*x) u(0) = u(1) = 0
    % The discretization of the problem is -u_(j - 1) + 
    % (2 + h^2*delta^2) * u_j - u_(j + 1) = h^2 * sin(n * pi * x).  The
    % resulting tridiagonal matrix has the following rows,
    %       [ ... 0 -1 (2 + h^2 * delta^2) - 1 0 ... ]
    A_x_output = x_input;
    x_len = N;

    % REDACTED: matrix multiplication code is available on request

end