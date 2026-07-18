% This is the project manager for the conjugate gradient project which has
% three parts involving different initial conditions.  The first involves 
% the base conjugate gradient method, the second compares conjugate 
% gradient with preconditioned conjugate gradient, and the last compares 
% conjugate gradient, conjugate residual, and preconditioned conjugate 
% gradient methods.  Different preconditioning methods will be used for the 
% second and third section.  Each section will involve populating the
% initial guess, calling the solver, and determining the total convergence
% at the end.  In addition, since the A matrix is sparce for each of the
% problems there will be custom matrix multiplication modules for each.
% Moving forward conjugate gradient = CG and preconditioned conjugate
% gradient = PCG.

%% Question 3 - CG computationally
% The first step is to initialize the b vector which is going to random
% integers from 1 to 100
b = randi([1 100],1000,1);
epsilon = 10^-10;
initial_guess = zeros(1000, 1);

% Uncomment the solver to run
% [u_final, num_iterations] = Conjugate_Gradient(epsilon, initial_guess, b, 3, 0, 0);


%% Question 4 - CG / PCG with Jacobi and SSOR
% For this part f(x)  = -6e^x + 15xe^x + 18x^2e^x + 3x^3e^x
% Our first goal is to create the grid based on h and then populate the b
% vector so we can use CG and PCG methods to solve it.  The number of
% interior points = h^-1 - 1.  The total grid size is going to be the
% number of interior points + 2 since we have to account for the boundary
% conditions
h = .1;
num_interior_points = h^(-1) - 1;
num_grid_points = num_interior_points + 2;

b = zeros(num_grid_points, 1);
u_exact = zeros(num_grid_points, 1);
initial_guess = zeros(num_grid_points, 1);
cur_x = 0;
e = exp(1);
for i = 2:(num_grid_points - 1)
    cur_x = cur_x + h;
    e_x = e^cur_x;
    b(i) = h^2 * (3 * e_x * (3 * cur_x + cur_x^2));
    u_exact(i) = 3 * e_x * (cur_x - cur_x^2);
end

% For the problem we are going to move the h^2 term to the right hand side
% b = -b * h^2;

% Uncomment the solver blocks to run
% % The first test is using conjugate gradient
% [u_final_4_1, num_iterations_4_1, error_4_1] = Conjugate_Gradient(epsilon, initial_guess, b, 4, 0, 0)
% 
% norm(u_final_4_1 - u_exact)

% % The next test involves using PCG with jacobi
% [u_final_4_2, num_iterations_4_2, error_4_2] = Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, 1, 4, num_grid_points)
% 
% norm(u_final_4_2 - u_exact)

% % The last test involves using PCG with SSOR
% [u_final_4_3, num_iterations_4_3, error_4_3] = Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, 2, 4, num_grid_points)
% 
% diff = norm(u_final_4_3 - u_exact)

%% Question 5 - CG / CR / PCG with ILU / point SSOR / line SSOR
% This is a 2 dimensional problem but due to the structure of gradient
% methods the values have to be stored in vectors of length
% num_grid_points^2 instead of as a matrix.  With this change the boundary
% points are throughout the vector so cases need to be built into the loops
% that populate the vectors.
epsilon = 10^-7;
h = .005;
num_interior_points = h^(-1) - 1;
num_grid_points = num_interior_points + 2;
vector_length = num_grid_points^2;

b = zeros(vector_length, 1);
u_exact = zeros(vector_length, 1);
p = 10;

cur_x = 0;
cur_y = 0;
current_vector_index = num_grid_points;
e = exp(1);
for i = 2:(num_grid_points - 1)
    cur_x = cur_x + h;
    cur_y = 0;
    for j = 1:num_grid_points
        current_vector_index = current_vector_index + 1;
        if j == 1
            b(current_vector_index) = 0;
            u_exact(current_vector_index) = 0;
        elseif j == num_grid_points
            b(current_vector_index) = 0;
            u_exact(current_vector_index) = 0;
        else
            cur_y = cur_y + h;
            term_1 = 3 * e^cur_x * e^cur_y;
            term_2 = p * (3 * cur_x + cur_x^2) * (cur_y - cur_y^2);
            term_3 = (3 * cur_y + cur_y^2) * (cur_x - cur_x^2);
            b(current_vector_index) = h^2 * term_1 * (term_2 + term_3);
            u_exact(current_vector_index) = 3 * e^cur_x * e^cur_y *  ...
                (cur_x - cur_x^2) * (cur_y - cur_y^2);
        end
    end
end
initial_guess = zeros(vector_length, 1);

% Uncomment the solver blocks to run
% % [u_final_5_1, num_iterations_5_1, error_5_1] = Conjugate_Gradient(epsilon, initial_guess, b, 5, p, num_grid_points)
% %
% % norm(u_final_5_1 - u_exact)
% % 
% [u_final_5_2, num_iterations_5_2, error_5_2] = Conjugate_Residual(epsilon, initial_guess, b, 5, p, num_grid_points)
% 
% norm(u_final_5_2 - u_exact)
% % 
% [u_final_5_3, num_iterations_5_3, error_5_3] = Conjugate_Residual(epsilon, initial_guess, b, 5, p, num_grid_points)
% 
% norm(u_final_5_3 - u_exact)
% % 
% [u_final_5_4, num_iterations_5_4, error_5_4] = Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, 3, 5, vector_length, p, num_grid_points)
% 
% norm(u_final_5_4 - u_exact)
% % 
% [u_final_5_5, num_iterations_5_5, error_5_5] = Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, 4, 5, vector_length, p, num_grid_points)
% 
% norm(u_final_5_5 - u_exact)
% %
% [u_final_5_6, num_iterations_5_6, error_5_6] = Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, 5, 5, vector_length, p, num_grid_points)
% 
% norm(u_final_5_6 - u_exact)