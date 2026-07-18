function [A_b] = Jacobi_Matrix_Vector_Multiplication(A,b)
    % This function multiplies a diagonal matrix with a vector
    % Inputs:
    %   A = matrix
    %   b = vector
    % Outputs:
    %   A_b = product of A and b
    A_b = zeros(size(A, 1));

    for i = 1:range(size(A, 1))
        A_b(i) = A(i, i) * b(i);
    end
end