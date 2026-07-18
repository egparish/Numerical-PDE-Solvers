function [Ax] = Question_5_Matrix_Multiplication(x, p, num_grid_points)
    % This function is going to handle matrix vector multiplications for
    % part 3, this is a 2 dimensional problem so it is going to be more
    % complicated than the prior two matrix multiplication functions.
    % From the structure of the problem we want to avoid all of the
    % boundary points and just do the matrix vector multiplication for the
    % interior points.  From the discretization each interior point has 5
    % terms that combine to form the dot product.  The rows have the
    % following form
    %            [0 ... -1 ... 0 -p 2(p + 1) -p 0 ... -1 ... 0]
    % where the -1 values are i +- num_grid_points where i is the current
    % row.  
    % Inputs:
    %   x = the vector we are multiplying with A
    %   p = parameter from the problem
    %   num_grid_points = if the vector was a matrix this is how long each
    %   row is, this includes the 2 boundary points
    % Outputs:
    %   Ax = the matrix vector product
    size_x = size(x, 1);
    Ax = zeros(size_x, 1);
    diag = 2 * (p + 1);

    % REDACTED: part 3 matrix multiplication code is available on request
end