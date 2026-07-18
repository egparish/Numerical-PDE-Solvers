function [Ax] = Question_3_Matrix_Multiplication(x)
    % This function is going to handle matrix vector multiplications for
    % part 1.  This problem involves a 1000x1000 diagonal matrix where
    % the diagonal values are 1, 2, ..., 9, 10 and they repeat.  This
    % function is going to utilize A being sparce
    % Inputs:
    %   x = the vector we are multiplying with A
    % Outputs:
    %   Ax = the matrix vector product
    Ax = x;
    diagonal_value = 1;
    for i = 1:length(x)
        Ax(i) = diagonal_value * x(i);

        % We want to update the diagonal value, if it gets to 11, we want
        % to set it to 1
        diagonal_value = diagonal_value + 1;
        if diagonal_value == 11
            diagonal_value = 1;
        end
    end
end