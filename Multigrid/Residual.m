function [residual_h] = Residual(b, x, delta, h, N)
    % This function is going to calculate the residual when moving down the
    % V-cycle.  The residual = b - A^h*X.  For this problem, A is a sparce
    % matrix whose structure depends on the discretization of the problem.
    % For efficient computing, the matrix multiplication function will be 
    % used.
    % Inputs:
    %   b = the answer to the boundary value problem
    %   x = the current x values
    %   delta = delta from boundary value problem
    %   h = grid spacing
    %   N = size of grid
    % Outputs
    %   residual_h = the current residual vector

    % The residual vector is the same size as the b and x vectors so we are
    % going to copy it and update all of the values
    residual_h = b;
    A_x = Matrix_Multiplication(x, delta, h, N);

    residual_len = N;
    for i = 1:residual_len
        residual_h(i) = b(i) - A_x(i);
    end
end