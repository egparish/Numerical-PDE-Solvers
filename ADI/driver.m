% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

% This is the driver for the ADI solver, it handles retrieving the data,
% running the ADI loop, and saving the data.  This will have code for both
% parts, with the desired code needed to be uncommented to run.

% j = x coordinate iteration
% k = y coordinate iteration

% The first step is to setup the scheme using the setup function which
% will read the needed parameters and generates the bands from the 2D
% stencil from ADI
% % For Part 1
% [NX, NY, XL, XR, YBOT, YTOP, ALPHA, BETA, D1, D2, DT, ...
%     TFIN, NSTEPS, DTPLOT, AX, BX, CX, AY, BY, CY] = setup();

% % For Part 2
% [NX, NY, XL, XR, YBOT, YTOP, SIGMA, Xc, Yc, D1, D2, DT, ...
%     TFIN, NSTEPS, DTPLOT, AX, BX, CX, AY, BY, CY] = setup();

% To write to the adi.dat file during the loop we need to open it
validate = fopen('adi.dat','w');

% Initialize the time variables
T = 0;
DTH = DT / 2;
TPLOT = 0;
h_x = 1 / (NX - 1);
h_y = 1 / (NY - 1);
iplot = 0;

% Initialize the solution grid with the initial conditions
% % Part 1
% u = init(NX, NY, ALPHA, BETA);

% % Part 2
% u = init(NX, NY, SIGMA, Xc, Yc);

% Initialize the right hand side vector that will be used with the Thomas
% algorithm
RHS = zeros(NX, 1);

% These are going to store the temporary solution
UHALF = u;
UNEW = u;

% This is the loop for the solver which is going to be comprised of two
% parts, the first being a traditional ADI solver and the second being two
% steps of backwards Euler when D2 = 0.  This will allow accuracy
% comparisons between ADI and the double backwards Euler.  This loop will
% also take care of saving the data periodically for the report.  For the
% ADI section it will be comprised of two steps, the first being an
% implicit step in x and an explicit step in y, with the second being an 
% implicit step in y and an explicit step in x.  Each will take care of 
% populating the right hand side using the discretization, enforcing the 
% boundary conditions, and then solving the system using a tridiagonal 
% solver.  The backward Euler step will involve populating the right hand
% side using the prior solution and then using the tridiagonal solver.

% The first loop cycles through the solution time
for ITIME = 1:NSTEPS
    
    % REDACTED: the solver loop is available on request

end  % END ITIME LOOP

% Close the adi.dat file after we are done
fclose(validate);

% The code to calculate the exact solution, save the data, and create plots
% has been removed