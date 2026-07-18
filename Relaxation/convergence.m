function [error] = convergence(x1,x2)
    % This function calculate the Manhattan or L1 norm between two vectors
    % Inputs:
    %   x1 = first vector
    %   x2 = second vector
    % Outputs:
    %   error = total distance

    error = 0;

    for i = 1:size(x1)
        error = error + abs(x1(i) - x2(i));
    end
end