% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u_final, num_iterations, error] = Conjugate_Gradient(epsilon, ...
    initial_guess, b, problem, p, num_grid_points)
    % This function is going to implement the conjugate gradient method to
    % solve a system of equations in the form Ax=b.  
    % Inputs:
    %   epsilon = convergence criteria
    %   initial_guess = the value of x_0
    %   b = solution to the system
    %   problem = the question conjugate gradient is being used for, it is
    %   going to be 3, 4, or 5.  It will determine the matrix
    %   multiplication function being used when calculating Ap_k
    %   p = parameter for question 5
    %   num_grid_points = size of the rows of the 2D vector, only used for 
    %   problem 5 will be 0 otherwise
    % Outputs:
    %   u_final = the final u values
    %   num_iterations = the number of iterations needed for the method to
    %   converge
    %   error = the final CG error

    % As mentioned in the project manager separate functions will be used
    % to handle the Ax calculations depending on the problem.  After this
    % the needed vectors will be initialized and then the solver will be run
    % until it converges.  Any vector ending in _old refer to the (k - 1)
    % vector while the vectors ending in _new refer to the (k + 1) vectors.
    % The value of r_initial does not change since it is used for the
    % convergence check.
    x_old = initial_guess;

    if problem == 4
        Ax_0 = Question_4_Matrix_Multiplication(x_old);
    else
        Ax_0 = Question_5_Matrix_Multiplication(x_old, p, num_grid_points);
    end

    % REDACTED: CG code is available on request

    % Update the return value with the result from the conjugate gradient
    % solver
    u_final = x_new;
end