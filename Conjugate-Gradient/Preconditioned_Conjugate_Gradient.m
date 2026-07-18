% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, num_iterations, error] = ...
    Preconditioned_Conjugate_Gradient(epsilon, initial_guess, b, ...
    conditioner, problem, N, p, num_grid_points)
    % This function is going to implement the preconditioned conjugate 
    % gradient method to solve a system of equations in the form Ax=b.  
    % Inputs:
    %   epsilon = convergence criteria
    %   initial_guess = the value of x_0
    %   b = solution to the system
    %   conditioner = boolean value (1) if we are using Jacobi as the
    %   preconditioner, if it is 0 then we are going to use SSOR
    %   problem = the question preconditioned conjugate gradient is being 
    %   used for, it is going to be 4 or 5.  It will determine the matrix
    %   multiplication function being used when calculating Ap_k
    %   N = number of equations including the boundary points
    %   p = parameter for question 5
    %   num_grid_points = size of the rows of the 2D vector, only used for 
    %   problem 5 will be 0 otherwise
    % Outputs:
    %   u_final = the final u values
    %   num_iterations = the number of iterations needed for the method to
    %   converge
    %   error = the error of the final PCG

    % As mentioned in the project manager separate functions will be used
    % to handle the Ax calculations depending on the problem.  After this
    % the needed vectors will be initialized and then the solver will be run
    % until it converges.  Any vector ending in _old refer to the (k - 1)
    % vector while the vectors ending in _new refer to the (k + 1) vectors.
    % The value of r_initial does not change since it is used for the
    % convergence check.
    x_old = initial_guess;
    omega = 1;

    r_initial = b;
    r_old = r_initial;
    x_new = x_old;
    r_new = r_old;

    % z_old = M^-1 r_0 - The preconditioner being used depends on the
    % values of the Jacobi variable, if it is 1 then we use Jacobi and
    % otherwise we are going to use SSOR
    if conditioner == 1
        % If we are using Jacobi, we are using the diagonals of A as M
        z_old = Jacobi_Preconditioner(r_old, N);
    elseif conditioner == 2
        % If we get here we are going to solve this used two SOR sweeps,
        % one of them is down and the other is up to solve Mz_0 = r_0
        z_old = SSOR_Preconditioner(omega, r_old, N);
    elseif conditioner == 3
        % If we get here this is for question 5 and we are going to use the
        % modified SSOR preconditioner
        z_old = Q5_Point_SSOR_Preconditioner(omega, r_old, N, p, num_grid_points);
    elseif conditioner == 4
        z_old = Q5_Line_SSOR_Preconditioner(omega, r_old, N, p, num_grid_points);
    elseif conditioner == 5
        z_old = Q5_ILU_Preconditioner(r_old, N, p, num_grid_points);
    end

    p_new = z_old;

    ZTR = Dot_Product(z_old, r_old);

    % We are going to run a while loop until the method converges or when
    % ||r_new|| / ||r_initial|| < epsilon
    r_initial_norm = norm(r_initial);
    error = 1;
    num_iterations = 0;

    % REDACTED: PCG code is available on request

    % Update the return value with the result from the conjugate gradient
    % solver
    u_final = x_new;
end