% This is going to be the project manager for project 2 which covers
% multigrids.  To split this code into different parts the following modules
% will be needed:
%   1) Relaxation Module for any Grid
%   2) Matrix Multiplication, specific to PDE
%   3) Residual
%   4) Fine -> Coarse Grid Transfer
%   5) Coarse -> Fine Grid Transfer
%   6) Coarse Grid Correction
%   7) Weighted Jacobi
%   8) Gaussian Elimination
% There will also be two types of data structures, one for V cycles and a
% separate one for W cycles

% Each grid has 2 boundary points and N interior points for a total of N +
% 2 points for the grid.  For this project the number of grids is going to
% be based off of the coarsest grid by saying that N + 1 = alpha * 2^Q
% where Q is the total number of grids minus 1 and alpha = N + 1 where N is
% the number of interior points of the smallest grid. For the h values at 
% each grid, h = 1 / (N + 1).  When moving down to coarser grids the odd 
% indexed points will be dropped.  The first step is to set the values and
% the initialize and fill the vectors that will be used.
alpha = 4;
num_grids = 6;
delta = 1;
n = 4;
epsilon = 10^-5;
Q = num_grids;
omega = 0.67;
max_grid_size = alpha * 2^(num_grids - 1) + 1;
h_start = 1 / (max_grid_size - 1);
total_cost = 0;
v_down = 10;
v_up = 1000;

% Assign the b values using the discretization of the problem and
% maintaining the Dirichlet boundary conditions.
b = zeros(1, max_grid_size);
first_term = n^2 * pi^2;
second_term = (sin(n * pi * h_start / 2))^2;
third_term = 1/4 * (n^2 * pi^2 * h_start^2);
C_tilda = 1 / ((first_term * second_term / third_term) + delta);
% first_term_2 = (sin(n * pi * h_start / 2)^2);
% second_term_2 = 4 / h_start^2;
% C_tilda = 1 / (first_term_2 * second_term_2 + delta);
%C_tilda = h_start^2 / (4 * sin(n * pi * h_start)^2 + delta * h_start^2);

for i = 1:(max_grid_size)
    if i == 1 
        b(i) = 0;
    elseif i == max_grid_size
        b(max_grid_size) = 0;
    else
        % Since we start at 0, x_j = 0 * i * h_start
        b(i) = C_tilda * sin(n * pi * i * h_start) * h_start^2;
        % b(i) = sin(n + pi * i * h_start) * h_start^2;
    end
end

% Once we have everything initialized we can solve the system using
% V-cycles which is implemented in a separate function.  The following
% parameters need to be passed to that function:
%   b
%   start_u
%   alpha
%   num_grids
%   delta
%   n
%   v_down
%   v_up
%   epsilon
% From here, all of the vectors and values are initialized and populated
% which allowed the V-cycle pass to be complete and be run until the entire
% solution converges to the desired tolerance or maximum cycles are run.
% There are different versions so uncomment the part needed
error = 1;
total_cost = 0;
tol = 10^-5;
start_u = zeros(1, max_grid_size);
iter = 0;

% while error > tol & iter < 10
%     iter = iter + 1
%     start_u;
%     [u_final, error_iter, computational_cost] = V_Cycle(b, start_u, alpha, ...
%         num_grids, delta, n, v_down, v_up, epsilon);
%     % [u_final, error_iter, computational_cost] = V_Cycle(start_u);
%     error = error_iter
%     start_u = u_final;
%     b;
%     total_cost = total_cost + computational_cost
% end

% while error > tol & iter < 3
%     iter = iter + 1
%     start_u;
% 
%     % % Choose between W or V cycles
%     % [u_final, error_iter, computational_cost] = W_Cycle(b, start_u, alpha, ...
%     %     num_grids, delta, n, v_down, v_up, epsilon);
%     % [u_final, error_iter, computational_cost] = V_Cycle(start_u);
% 
%     error = error_iter
%     start_u = u_final;
%     b;
%     total_cost = total_cost + computational_cost
% end
