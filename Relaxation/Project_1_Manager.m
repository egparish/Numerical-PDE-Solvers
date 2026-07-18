% We are considering the discretization of the following boundary value
% problem: -(epsilon * u_xx + u_yy) = 0 for 0 <= x <= 1, 0 <= y <= 1 a,
% with u = 0 on all boundaries.  We are going to discretize as described in
% the notes for the following problems:
%   a) u_j,k = (-1)^(jk)
%   b) u_j,k = sin(8πjh)sin(8πkh)
%   c) u_j,k = sin(πjh) sin(πkh)
% This code to solve these problems has been modularized so each sub part
% of the program is in its own file.  This file is going to manage the
% project.  It is going to generate the A matrix and b vectors along with
% calling the various relaxation methods and handling the plots that we
% need to generate.

% The first step is to generate the A matrix, due to the methods we are
% using we want to include the boundary points as a part of the A matrix.
% If we are using a specified N value then the resulting A matrix is going
% to be (N+1) x (N+1).  These added values are only going to be used for
% the discretization of the problem.  This is a 2D problem and as a result
% the discretization being used will be 

% N is used to determine the number of interior points on the grid
N = 1;

% For the discretization of the PDE I am going to use the following:
%   u_xx * H^2 = u_(j - 1, k) - 2_(j, k) + u_(j + 1, k)
%   u_yy * H^2 = u_(j, k - 1) - 2_(j, k) + u_(, k + 1)
% That leaves us with
%   0 = -(epsilon(u_(j - 1, k) - 2_(j, k) + u_(j + 1, k)) - 
%                 u_(j, k - 1) - 2_(j, k) + u_(, k + 1))



% This is a Dirichlet problem, all of the boundary conditions are 0

% The discretization, u_(j, k) and eplison give us the final A.  The
% discretization gives us the form of A while u_(j, k) and eplison scales
% it

% From the discretization and boundary conditions of the problem we know
% that the base shape of A is going to be ...

% the iths row will be [... 1   0   0   1   -4   1   0   0   1 ...]
%                          i-4         i-1   i  i+1         i+4


% A is going to be similar to the matrix on page 18 of the notes but since
% all we care about is Ax then we can just compute a modified matrix
% multiplecation

% For this program, we are going to store A as a sparce matrix?  N is small
% enough in this case but in other situations I would store it differently.
%  The matrix vector multiplecation Ax is going to be computed efficiently
%  since A is sparce and will be done from a separate function.  A will be
%  (N+2) by (N+2) to account for the boundary points, in this case they are
%  all 0 but it is easier for me to think about it this way.

% I am going to test the programs on a system that I know works.  In this
% case the system is going to be 

A_test = [0 1 1 1;
          1 0 1 1;
          1 1 0 1;
          1 1 1 0];

b_test = [70; 75; 80; 75];

test_ans = [30; 25; 20; 25];

% 
% % This is test code for the Jacobi Matrix Multiplier, I am going to do a
% % 10*10 matrix * Eye(10) * 2
% A_mat_test = [-4, 1, 0, 0, 1, 0, 0, 0, 0, 0;
%               1, -4, 1, 0, 0, 1, 0, 0, 0, 0;
%               0, 1, -4, 1, 0, 0, 1, 0, 0, 0;
%               0, 0, 1, -4, 1, 0, 0, 1, 0, 0;
%               1, 0, 0, 1, -4, 1, 0, 0, 1, 0;
%               0, 1, 0, 0, 1, -4, 1, 0, 0, 1;
%               0, 0, 1, 0, 0, 1, -4, 1, 0, 0;
%               0, 0, 0, 1, 0, 0, 1, -4, 1, 0;
%               0, 0, 0, 0, 1, 0, 0, 1, -4, 1;
%               0, 0, 0, 0, 0, 1, 0, 0, 1, -4];
% 
% D_mat_test = [.5, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%               0, .5, 0, 0, 0, 0, 0, 0, 0, 0;
%               0, 0, .5, 0, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, .5, 0, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, .5, 0, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, .5, 0, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, .5, 0, 0, 0;
%               0, 0, 0, 0, 0, 0, 0, .5, 0, 0;
%               0, 0, 0, 0, 0, 0, 0, 0, .5, 0;
%               0, 0, 0, 0, 0, 0, 0, 0, 0, .5];
% 
% true_val = D_mat_test * A_mat_test;
% 
% my_val = Jacobi_Matrix_Multiplication(A_mat_test, D_mat_test);
% 
% true_val - my_val
% 
% % This is going to test how I take the inverse of a diagonal matrix
% diag_mat_test = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%                  0, 2, 0, 0, 0, 0, 0, 0, 0, 0;
%                  0, 0, 3, 0, 0, 0, 0, 0, 0, 0;
%                  0, 0, 0, 4, 0, 0, 0, 0, 0, 0;
%                  0, 0, 0, 0, 5, 0, 0, 0, 0, 0;
%                  0, 0, 0, 0, 0, 6, 0, 0, 0, 0;
%                  0, 0, 0, 0, 0, 0, 7, 0, 0, 0;
%                  0, 0, 0, 0, 0, 0, 0, 8, 0, 0;
%                  0, 0, 0, 0, 0, 0, 0, 0, 9, 0;
%                  0, 0, 0, 0, 0, 0, 0, 0, 0, 10];
% 
% diag_inv_test = inv(diag_mat_test)
% 
% % This is going to test the Jacobi_G_Matrix_X_Vector_Multiplecation(G,X)
% X_vect_test = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
% 
% funct_test = Jacobi_G_Matrix_X_Vector_Multiplecation(A_mat_test, X_vect_test)
% 
% A_mat_test * X_vect_test

% I am just stupid.  We are solving for the u's in the u_j,k initial data
% using the discretization described.  We are not storing any part of A, we
% are constantly using the new values of u in the equations.

% The solutions at each step are the u values, we store only store the u
% values and use the discretization matrix.  We are including boundary
% points


% This is going to be code to generate the initial data and test the
% algorithms
N = 16;
h = 1 / (N + 1);
epsilon = 1
u_matrix_a_16 = zeros(N + 2, N + 2);
u_matrix_b_16 = zeros(N + 2, N + 2);
u_matrix_c_16 = zeros(N + 2, N + 2);

% Populate the matrix
for j = 2:(N + 1)
    for k = 2:(N + 1)
        u_matrix_a_16(j, k) = (-1)^(j + k);

        u_matrix_b_16(j, k) = sin(8 * pi * j * h) * sin(8 * pi * k * h);

        u_matrix_c_16(j, k) = sin(pi * j * h) * sin(pi * k * h);
    end
end

N = 32;
h = 1  / (N + 1);
u_matrix_a_32 = zeros(N + 2, N + 2);
u_matrix_b_32 = zeros(N + 2, N + 2);
u_matrix_c_32 = zeros(N + 2, N + 2);
for j = 2:(N + 1)
    for k = 2:(N + 1)
        u_matrix_a_32(j, k) = (-1)^(j + k);

        u_matrix_b_32(j, k) = sin(8 * pi * j * h) * sin(8 * pi * k * h);

        u_matrix_c_32(j, k) = sin(pi * j * h) * sin(pi * k * h);
    end
end


N = 64;
h = 1  / (N + 1);
u_matrix_a_64 = zeros(N + 2, N + 2);
u_matrix_b_64 = zeros(N + 2, N + 2);
u_matrix_c_64 = zeros(N + 2, N + 2);
for j = 2:(N + 1)
    for k = 2:(N + 1)
        u_matrix_a_64(j, k) = (-1)^(j + k);

        u_matrix_b_64(j, k) = sin(8 * pi * j * h) * sin(8 * pi * k * h);

        u_matrix_c_64(j, k) = sin(pi * j * h) * sin(pi * k * h);
    end
end


N = 128;
h = 1  / (N + 1);
u_matrix_a_128 = zeros(N + 2, N + 2);
u_matrix_b_128 = zeros(N + 2, N + 2);
u_matrix_c_128 = zeros(N + 2, N + 2);
for j = 2:(N + 1)
    for k = 2:(N + 1)
        u_matrix_a_128(j, k) = (-1)^(j + k);

        u_matrix_b_128(j, k) = sin(8 * pi * j * h) * sin(8 * pi * k * h);

        u_matrix_c_128(j, k) = sin(pi * j * h) * sin(pi * k * h);
    end
end

% Uncomment the section to run
% [u_final_a_16, iterations_a_16, max_val_a_16] = Point_Jacobi(u_matrix_a_16, epsilon, 16);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Point_Jacobi(u_matrix_b_16, epsilon, 16);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Point_Jacobi(u_matrix_c_16, epsilon, 16);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Point_Jacobi(u_matrix_a_32, epsilon, 32);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Point_Jacobi(u_matrix_b_32, epsilon, 32);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Point_Jacobi(u_matrix_c_32, epsilon, 32);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Point_Jacobi(u_matrix_a_64, epsilon, 64);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Point_Jacobi(u_matrix_b_64, epsilon, 64);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Point_Jacobi(u_matrix_c_64, epsilon, 64);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Point_Jacobi(u_matrix_a_128, epsilon, 128);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Point_Jacobi(u_matrix_b_128, epsilon, 128);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Point_Jacobi(u_matrix_c_128, epsilon, 128);

% [u_final_a_16, iterations_a_16, max_val_a_16] = Weighted_Jacobi(u_matrix_a_16, epsilon, 16, 0.9);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Weighted_Jacobi(u_matrix_b_16, epsilon, 16, 0.9);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Weighted_Jacobi(u_matrix_c_16, epsilon, 16, 0.9);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Weighted_Jacobi(u_matrix_a_32, epsilon, 32, 0.9);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Weighted_Jacobi(u_matrix_b_32, epsilon, 32, 0.9);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Weighted_Jacobi(u_matrix_c_32, epsilon, 32, 0.9);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Weighted_Jacobi(u_matrix_a_64, epsilon, 64, 0.9);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Weighted_Jacobi(u_matrix_b_64, epsilon, 64, 0.9);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Weighted_Jacobi(u_matrix_c_64, epsilon, 64, 0.9);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Weighted_Jacobi(u_matrix_a_128, epsilon, 128, 0.9);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Weighted_Jacobi(u_matrix_b_128, epsilon, 128, 0.9);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Weighted_Jacobi(u_matrix_c_128, epsilon, 128, 0.9);


% [u_final_a_16, iterations_a_16, max_val_a_16] = Gauss_Seidel(u_matrix_a_16, epsilon, 16);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Gauss_Seidel(u_matrix_b_16, epsilon, 16);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Gauss_Seidel(u_matrix_c_16, epsilon, 16);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Gauss_Seidel(u_matrix_a_32, epsilon, 32);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Gauss_Seidel(u_matrix_b_32, epsilon, 32);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Gauss_Seidel(u_matrix_c_32, epsilon, 32);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Gauss_Seidel(u_matrix_a_64, epsilon, 64);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Gauss_Seidel(u_matrix_b_64, epsilon, 64);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Gauss_Seidel(u_matrix_c_64, epsilon, 64);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Gauss_Seidel(u_matrix_a_128, epsilon, 128);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Gauss_Seidel(u_matrix_b_128, epsilon, 128);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Gauss_Seidel(u_matrix_c_128, epsilon, 128);

% [u_final_a_16, iterations_a_16, max_val_a_16] = Red_Black_Gauss_Seidel(u_matrix_a_16, epsilon, 16);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Red_Black_Gauss_Seidel(u_matrix_b_16, epsilon, 16);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Red_Black_Gauss_Seidel(u_matrix_c_16, epsilon, 16);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Red_Black_Gauss_Seidel(u_matrix_a_32, epsilon, 32);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Red_Black_Gauss_Seidel(u_matrix_b_32, epsilon, 32);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Red_Black_Gauss_Seidel(u_matrix_c_32, epsilon, 32);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Red_Black_Gauss_Seidel(u_matrix_a_64, epsilon, 64);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Red_Black_Gauss_Seidel(u_matrix_b_64, epsilon, 64);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Red_Black_Gauss_Seidel(u_matrix_c_64, epsilon, 64);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Red_Black_Gauss_Seidel(u_matrix_a_128, epsilon, 128);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Red_Black_Gauss_Seidel(u_matrix_b_128, epsilon, 128);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Red_Black_Gauss_Seidel(u_matrix_c_128, epsilon, 128);

% omega_SSOR = 1.9;
% omega_LINE_SSOR = 1.5;
% epsilon = .001;
% [u_comparisson_SSOR, iterations_SSOR, max_val_SSOR] = SSOR(u_matrix_c_128, epsilon, 128, omega_SSOR);
% [u_comparisson_LINE_SSOR, iterations_LINE_SSOR, max_val_LINE_SSOR] = Line_SSOR(u_matrix_c_128, epsilon, 128, omega_LINE_SSOR);


% [u_final_a_16, iterations_a_16, max_val_a_16] = SOR(u_matrix_a_16, epsilon, 16, omega);
% [u_final_b_16, iterations_b_16, max_val_b_16] = SOR(u_matrix_b_16, epsilon, 16, omega);
% [u_final_c_16, iterations_c_16, max_val_c_16] = SOR(u_matrix_c_16, epsilon, 16, omega);
% [u_final_a_32, iterations_a_32, max_val_a_32] = SOR(u_matrix_a_32, epsilon, 32, omega);
% [u_final_b_32, iterations_b_32, max_val_b_32] = SOR(u_matrix_b_32, epsilon, 32, omega);
% [u_final_c_32, iterations_c_32, max_val_c_32] = SOR(u_matrix_c_32, epsilon, 32, omega);
% [u_final_a_64, iterations_a_64, max_val_a_64] = SOR(u_matrix_a_64, epsilon, 64, omega);
% [u_final_b_64, iterations_b_64, max_val_b_64] = SOR(u_matrix_b_64, epsilon, 64, omega);
% [u_final_c_64, iterations_c_64, max_val_c_64] = SOR(u_matrix_c_64, epsilon, 64, omega);
% [u_final_a_128, iterations_a_128, max_val_a_128] = SOR(u_matrix_a_128, epsilon, 128, omega);
% [u_final_b_128, iterations_b_128, max_val_b_128] = SOR(u_matrix_b_128, epsilon, 128, omega);
% [u_final_c_128, iterations_c_128, max_val_c_128] = SOR(u_matrix_c_128, epsilon, 128, omega);

% [u_final_a_16, iterations_a_16, max_val_a_16] = SSOR(u_matrix_a_16, epsilon, 16, omega);
% [u_final_b_16, iterations_b_16, max_val_b_16] = SSOR(u_matrix_b_16, epsilon, 16, omega);
% [u_final_c_16, iterations_c_16, max_val_c_16] = SSOR(u_matrix_c_16, epsilon, 16, omega);
% [u_final_a_32, iterations_a_32, max_val_a_32] = SSOR(u_matrix_a_32, epsilon, 32, omega);
% [u_final_b_32, iterations_b_32, max_val_b_32] = SSOR(u_matrix_b_32, epsilon, 32, omega);
% [u_final_c_32, iterations_c_32, max_val_c_32] = SSOR(u_matrix_c_32, epsilon, 32, omega);
% [u_final_a_64, iterations_a_64, max_val_a_64] = SSOR(u_matrix_a_64, epsilon, 64, omega);
% [u_final_b_64, iterations_b_64, max_val_b_64] = SSOR(u_matrix_b_64, epsilon, 64, omega);
% [u_final_c_64, iterations_c_64, max_val_c_64] = SSOR(u_matrix_c_64, epsilon, 64, omega);
% [u_final_a_128, iterations_a_128, max_val_a_128] = SSOR(u_matrix_a_128, epsilon, 128, omega);
% [u_final_b_128, iterations_b_128, max_val_b_128] = SSOR(u_matrix_b_128, epsilon, 128, omega);
% [u_final_c_128, iterations_c_128, max_val_c_128] = SSOR(u_matrix_c_128, epsilon, 128, omega);

% [u_final_a_16, iterations_a_16, max_val_a_16] = Line_SSOR(u_matrix_a_16, epsilon, 16, 1.8);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Line_SSOR(u_matrix_b_16, epsilon, 16, 1.8);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Line_SSOR(u_matrix_c_16, epsilon, 16, 1.8);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Line_SSOR(u_matrix_a_32, epsilon, 32, 1.8);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Line_SSOR(u_matrix_b_32, epsilon, 32, 1.8);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Line_SSOR(u_matrix_c_32, epsilon, 32, 1.8);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Line_SSOR(u_matrix_a_64, epsilon, 64, 1.8);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Line_SSOR(u_matrix_b_64, epsilon, 64, 1.8);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Line_SSOR(u_matrix_c_64, epsilon, 64, 1.8);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Line_SSOR(u_matrix_a_128, epsilon, 128, 1.8);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Line_SSOR(u_matrix_b_128, epsilon, 128, 1.8);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Line_SSOR(u_matrix_c_128, epsilon, 128, 1.8);

% [u_final_a_16, iterations_a_16, max_val_a_16] = Kaczmarz(u_matrix_a_16, epsilon, 16, 25000);
% [u_final_b_16, iterations_b_16, max_val_b_16] = Kaczmarz(u_matrix_b_16, epsilon, 16, 25000);
% [u_final_c_16, iterations_c_16, max_val_c_16] = Kaczmarz(u_matrix_c_16, epsilon, 16, 25000);
% [u_final_a_32, iterations_a_32, max_val_a_32] = Kaczmarz(u_matrix_a_32, epsilon, 32, 25000);
% [u_final_b_32, iterations_b_32, max_val_b_32] = Kaczmarz(u_matrix_b_32, epsilon, 32, 25000);
% [u_final_c_32, iterations_c_32, max_val_c_32] = Kaczmarz(u_matrix_c_32, epsilon, 32, 25000);
% [u_final_a_64, iterations_a_64, max_val_a_64] = Kaczmarz(u_matrix_a_64, epsilon, 64, 25000);
% [u_final_b_64, iterations_b_64, max_val_b_64] = Kaczmarz(u_matrix_b_64, epsilon, 64, 25000);
% [u_final_c_64, iterations_c_64, max_val_c_64] = Kaczmarz(u_matrix_c_64, epsilon, 64, 25000);
% [u_final_a_128, iterations_a_128, max_val_a_128] = Kaczmarz(u_matrix_a_128, epsilon, 128, 25000);
% [u_final_b_128, iterations_b_128, max_val_b_128] = Kaczmarz(u_matrix_b_128, epsilon, 128, 25000);
% [u_final_c_128, iterations_c_128, max_val_c_128] = Kaczmarz(u_matrix_c_128, epsilon, 128, 25000);

% u_matrix_1;
% u_matrix_2;
% u_matrix_3;
% 
% max(max(u_matrix_3))

% [u_final_1, iterations_1] = Point_Jacobi(u_matrix_1, epsilon, N)

% [u_final, iterations_2] = Weighted_Jacobi(u_matrix, epsilon, N, .8)

% [u_final_3, iterations_3] = Gauss_Seidel(u_matrix, epsilon, N)

% [u_final_4, iterations_4] = Red_Black_Gauss_Seidel(u_matrix, epsilon, N)

% [u_final_5, iterations_5] = SOR(u_matrix, epsilon, N, 1.8)

% [u_final_6, iterations_6] = SSOR(u_matrix, epsilon, N, 1.8)

% [u_final_7, iterations_7] = Line_SSOR(u_matrix, epsilon, N, 1.8)

% [u_final_8, iterations_8] = Kaczmarz(u_matrix_3, epsilon, N, 25000)