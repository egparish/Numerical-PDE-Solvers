% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [x_final] = Thomas_Algorithm(epsilon, N, u, row)
    % This function is going to implement the Thomas algorithm for a
    % tridiagonal matrix, which involves a forward sweep to eliminate the
    % subdiagonal and a second sweep to solve the equation using back
    % substitution.  The A matrix is defined using epsilon the 
    % discretization of the original problem and the b values will be
    % derived from the u matrix
    % Inputs:
    %   epsilon = from the original PDE
    %   N = the number of data points in u, does not include the imaginary
    %   points, this is the number of equations we are solving
    %   u = the u values for the system
    %   row = row we are working with, used for u
    % Outputs
    %   x_final = the results of the relaxation method

    % For all rows of A the diagonal is going to be (2 + 2 * epsilon) and
    % the sub and super diagonal are -1.  It is important to remember when
    % dealing with u is that it is (N+2) x (N+2) so all indices used need 1
    % added to them.

    % For the ith element in the jth row we are solving for, b_i:
    %   b_i = u_(i, j + 1) + u(i, j - 1)

    % Define the ouput and scratch vector
    x_final = zeros(N, 1);
    scratch = zeros(N, 1);

    b = zeros(N, 1);
    for i = 1:N
        b(i) = u(i + 1, row + 1) + u(i + 1, row - 1);
    end
    
    % Coefficients
    a = -1; 
    c = -1; 
    d = 2 + 2 * epsilon;

    % REDACTED: Thomas algorithm code is available on request
end