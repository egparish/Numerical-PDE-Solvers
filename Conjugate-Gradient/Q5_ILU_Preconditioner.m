% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [z] = Q5_ILU_Preconditioner(b, N, p, num_grid_points)
    % This function is going to implement the Meijerink and Van Vorst 
    % variation of ILU decomposition to solve for the preconditioning
    % matrix.  M = LU and using the PCG algorithm this can be used to solve
    % 2 systems of equations, Ly_k = r_k and Uz_k = y_k from the original
    % equation, Mz_k = r_k
    % Inputs:
    %   b = solution to the system, includes boundary points
    %   N = number of equations, includes the boundary points
    %   p = parameter from original equation
    %   num_grid_points = equations per row
    % Outputs:
    %   z = the value of M^-1b

    % From the discretization of 4 we have the following for the rows of A
    %            [0 ... -1 ... 0 -p 2(p + 1) -p 0 ... -1 ... 0]
    % We want to perform an incomplete LU decomposition (ILU) on this matrix
    % using the Meijerink and Van Vorst variation which is going to
    % introduce extra terms but we are going to make the assumption for
    % easier computation.  We want to achieve the following form of L and U:
    %          | d_1  0   0   0   0  |
    %          | c_2 d_2  0   0   0  |
    %      L = |  0  c_3 d_3  0   0  |
    %          |  0   0  c_4 d_4  0  |
    %          | b_5  0   0  c_5 d_5 |
    %
    %          |  1  e_2  0   0  f_5 |
    %          |  0   1  e_3  0   0  |
    %      U = |  0   0   1  e_4  0  |
    %          |  0   0   0   1  e_5 |
    %          |  0   0   0   0   1  |
    %
    % Assuming A has the following form:
    %          | E_1 D_2  0   0  B_5 |
    %          | D_2 E_2 D_3  0   0  |
    %      A = |  0  D_3 E_3 D_4  0  |
    %          |  0   0  D_4 E_4 D_5 |
    %          | B_5  0   0  D_5 E_5 |
    %
    % We can solve the following system of equations, uppercase variables
    % are knows from the structure of A and the lowercase variables are in
    % the L and U matrices
    %   b_i = B_i                            i = n + 1, .., N
    %   c_i = D_i                            i = 2, ..., N
    %   d_i = E_i - D_i * e_i - B_i * f_i    i = 1, ..., N
    %   e_(i + 1) = D_(i + 1) / d_i          i = 1, ..., N
    %       or e_i = D_(i) / d_(i - 1)
    %   f_(i + n) = B_(i + n) / d_i          i = 1, ..., N

    % From here we have all of the information we need to solve for the
    % the L and U matrices.  Just the non zero values are going to be
    % stored since the future values depend on some of the old values.  As 
    % with the other preconditioners the boundary values are going to stay 
    % 0 and those rows, cols are going to be skipped but the structure of the
    % matrices will stay the same.  The first step is to solve for y using
    % b and L. 
    E = 2 * (p + 1);
    b_i = zeros(N, 1);
    c_i = zeros(N, 1);
    d_i = zeros(N, 1);
    e_i = zeros(N, 1);
    f_i = zeros(N, 1);

    % REDACTED: the ILU preconditioner code is available on request

end