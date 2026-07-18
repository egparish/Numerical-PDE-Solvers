function [A_B] = Jacobi_Matrix_Vector_Subtraction(A,B)
    % This function subtracts two diagonal matrices
    % Inputs:
    %   A = matrix
    %   B = matrix
    % Outputs:
    %   A_B = A - B
    A_B = eye(size(A), size(A));

    for i = 1:size(A)
        A_B(i, i) = A(i, i) - B(i, i);
    end
end