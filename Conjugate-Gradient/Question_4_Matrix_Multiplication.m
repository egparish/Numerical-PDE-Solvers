function [Ax] = Question_4_Matrix_Multiplication(x)
    % This function is going to handle matrix vector multiplications for
    % part 2.  This problem involves the discretization of -u_xx = f(x).
    %  The rows have the following form, 
    %                   [0 ... 0 -1 2 -1 0 ... 0]
    % with the first and last rows not being complete but these represent
    % the boundary points.  These two rows are going to be treated
    % differently but all of the interior calculations are the same.  
    % Inputs:
    %   x = the vector we are multiplying with A
    % Outputs:
    %   Ax = the matrix vector product
    Ax = zeros(size(x));

    % Handle the two edge cases outside of the loop
    Ax(1) = 0;
    Ax(end) = 0;

    for i = 2:(length(x) - 1)
        Ax(i) = -x(i - 1) + 2 * x(i) - x(i + 1);
    end
end