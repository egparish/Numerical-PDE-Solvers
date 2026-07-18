function [A_B] = Jacobi_Matrix_Matrix_Subtraction(I, A)
    % This function is going to compute I - D^-1A.  Since I only has values
    % on the diagonal, these are the only values that are important and the
    % rest can be built off of -A
    % Inputs:
    %   I = matrix
    %   A = matrix
    % Outputs:
    %   A_B = I - D^-1A

    A_B = -A;

    for i = 1:size(A, 1)
        A_B(i, i) = AI(i, i) + A(i, i);
    end
end