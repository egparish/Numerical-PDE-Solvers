% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

% This is setup function that will read the needed data from adiset, save
% it, write it to another file for validation reasons, and return the
% values to the file that calls the function.  This will have two versions
% of the function for each part with the desired version needing to be 
% uncommented to run.

% % Setup function for the first part
% function [NX, NY, XL, XR, YBOT, YTOP, ALPHA, BETA, D1, D2, DT, ...
%     TFIN, NSTEPS, DTPLOT, AX, BX, CX, AY, BY, CY] = setup()
%     % Inputs:
%     %   adiset (this is a stored data file)
%     % Outputs:
%     %   NX = number of points in the x direction
%     %   NY = number of points in the y direction
%     %   XL = left endpoint of the square
%     %   XR = right endpoint of the square
%     %   YBOT = bottom endpoint of the square
%     %   YTOP = top endpoint of the square
%     %   ALPHA = integer mode number for the IC
%     %   BETA = integer mode number for the IC
%     %   D1 = first diffusivities for the x direction
%     %   D2 = second diffusivities for the y direction
%     %   DT = timestep
%     %   TFIN = final time of computation
%     %   NSTEPS = maximum number of steps
%     %   DTPLOT = time intervals where the computed and exact solutions are
%     %            plotted
%     %   AX = lower diagonal for the x implicit step
%     %   BX = diagonal for the x implicit step
%     %   CX = upper diagonal for the x implicit step
%     %   AY = lower diagonal for the y implicit step
%     %   BY = diagonal for the y implicit step
%     %   CY = upper diagonal for the y implicit step
% 
%     % REDACTED: setup code is available on request
% end

% Setup function for the second part
% function [NX, NY, XL, XR, YBOT, YTOP, SIGMA, Xc, Yc, D1, D2, DT, ...
%     TFIN, NSTEPS, DTPLOT, AX, BX, CX, AY, BY, CY] = setup()
%     % Inputs:
%     %   adiset (this is a stored data file)
%     % Outputs:
%     %   NX = number of points in the x direction
%     %   NY = number of points in the y direction
%     %   XL = left endpoint of the square
%     %   XR = right endpoint of the square
%     %   YBOT = bottom endpoint of the square
%     %   YTOP = top endpoint of the square
%     %   SIGMA = problem parameter
%     %   Xc = X coordinate of the impulsive heating
%     %   Yc = Y coordinate of the impulsive heating
%     %   D1 = first diffusivities for the x direction
%     %   D2 = second diffusivities for the y direction
%     %   DT = timestep
%     %   TFIN = final time of computation
%     %   NSTEPS = maximum number of steps
%     %   DTPLOT = time intervals where the computed and exact solutions are
%     %   plotted
%     %   AX = lower diagonal for the x implicit step
%     %   BX = diagonal for the x implicit step
%     %   CX = upper diagonal for the x implicit step
%     %   AY = lower diagonal for the y implicit step
%     %   BY = diagonal for the y implicit step
%     %   CY = upper diagonal for the y implicit step
% 
%     % REDACTED: setup code is available on request
% end