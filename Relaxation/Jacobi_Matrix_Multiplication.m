function [DA] = Jacobi_Matrix_Multiplication(A,D)
    % This function computes the matrix product of a diagonal matrix D and
    % a matrix A.  The only computation is along the diagonal entries
    % Inputs:
    %   A = matrix
    %   D = diagonal matrix
    % Outputs:
    %   DA = product of D and A

    DA = eye(size(A));

    for i = 1:size(A)
        DA(i, i) = D(i, i) * A(i, i);
    end
end