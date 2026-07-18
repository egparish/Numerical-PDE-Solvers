% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [u] = macfb(u, n, lambda, c)
    % This function is going to complete one step of the MacCormack method
    % using a forwards predictor and a backwards corrector.
    % Inputs:
    %   u = current solution
    %   n = the granularity of the solution
    %   lambda = parameter for the solver
    %   c = function from the original problem
    % Outputs:
    %   u = the updated solution (u^(n + 1))

    % My preference is to preallocate the various vectors needed in the
    % function, f(j), ucap(j), fcap(j)
    f = zeros(1, n);
    ucap = zeros(1, n);
    fcap = zeros(1, n);

    % REDACTED: the 2-2 forwards-backward step is available on request
    
end
