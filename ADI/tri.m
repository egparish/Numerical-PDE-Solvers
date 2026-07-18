% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [SOL] = tri(A, B, C, D)
    % This function is going to implement the Thomas algorithm for a
    % tridiagonal matrix.  This will be used to solve for x from Ax = b
    % problems to avoid using the '\' operator and reduce the number of
    % operations needed.  The Thomas algorithm has two passes, the first to
    % eliminate the sup diagonal and the second to solve using back
    % substitution. For this implementation the inputs are going to be 
    % vectors with no matrices used.  To ensure the initial values are not 
    % changed there will be local vectors to hold the temp values.

    % Inputs:
    %   A = lower diagonal (length = n - 1)
    %   B = main diagonal (length = n)
    %   C = upper diagonal (length = n - 1)
    %   D = the right hand side (length n)
    % Outputs:
    %   SOL = the solution, x (length n)

    % REDACTED: the initialization code is available on request

end