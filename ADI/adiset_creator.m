% This script creates the adiset file used by the setup function.  There
% are two parts for this project, the first having the following initial
% conditions, u(0, x, y) = sin(alpha * pi * x) * sin(peta * pi * y), and
% the second having a heat pulse with Dirichlet boundary conditions.  Each
% part has their own parameter list with the desired part needing to be
% uncommented to run.

% % Parameters for part 1
% %       NX = number of points in the x direction
% %       NY = number of points in the y direction
% %       XL = left endpoint of the square
% %       XR = right endpoint of the square
% %       YBOT = bottom endpoint of the square
% %       YTOP = top endpoint of the square
% %       ALPHA = integer mode number for the IC
% %       BETA = integer mode number for the IC
% %       D1 = first diffusivities for the x direction
% %       D2 = second diffusivities for the y direction
% %       DT = timestep
% %       TFIN = final time of computation
% %       NSTEPS = maximum number of steps
% %       DTPLOT = time intervals where the computed and exact solutions are
% %       plotted
% 
% % Set the values
% NX = 51;
% NY = 51;
% 
% XL = 0; 
% XR = 1;
% YBOT = 0; 
% YTOP = 1;
% 
% ALPHA = 1;
% BETA = 1;
% 
% D1 = 1; 
% D2 = 0;
% 
% DT = 0.01;
% TFIN = 0.1;
% NSTEPS = 10000;
% DTPLOT = 0.05;
% 
% % --- Write the adiset file ---
% fid = fopen('adiset','w');
% 
% fprintf(fid,'%d %d\n', NX, NY);
% fprintf(fid,'%f %f\n', XL, XR);
% fprintf(fid,'%f %f\n', YBOT, YTOP);
% fprintf(fid,'%d %d\n', ALPHA, BETA);
% fprintf(fid,'%f %f\n', D1, D2);
% fprintf(fid,'%f\n', DT);
% fprintf(fid,'%f\n', TFIN);
% fprintf(fid,'%d\n', NSTEPS);
% fprintf(fid,'%f\n', DTPLOT);
% 
% fclose(fid);
% 
% disp('Created adiset file.');

% % Parameters for part 2
% %       NX = number of points in the x direction
% %       NY = number of points in the y direction
% %       XL = left endpoint of the square
% %       XR = right endpoint of the square
% %       YBOT = bottom endpoint of the square
% %       YTOP = top endpoint of the square
% %       SIGMA = problem parameter
% %       Xc = X coordinate of the impulsive heating
% %       Yc = Y coordinate of the impulsive heating
% %       D1 = first diffusivities for the x direction
% %       D2 = second diffusivities for the y direction
% %       DT = timestep
% %       TFIN = final time of computation
% %       NSTEPS = maximum number of steps
% %       DTPLOT = time intervals where the computed and exact solutions are
% %       plotted

% % Set the values
% NX = 4001;
% NY = 4001;
% 
% XL = 0; 
% XR = 1;
% YBOT = 0; 
% YTOP = 1;
% 
% SIGMA = 0.02;
% 
% Xc = 0.5;
% Yc = 0.5;
% 
% D1 = 1; 
% D2 = 1;
% 
% DT = 0.00005;
% TFIN = 0.01;
% NSTEPS = 100000000;
% DTPLOT = 0.001;
% 
% % --- Write the adiset file ---
% fid = fopen('adiset','w');
% 
% fprintf(fid,'%d %d\n', NX, NY);
% fprintf(fid,'%f %f\n', XL, XR);
% fprintf(fid,'%f %f\n', YBOT, YTOP);
% fprintf(fid,'%f\n', SIGMA);
% fprintf(fid,'%f %f\n', Xc, Yc);
% fprintf(fid,'%f %f\n', D1, D2);
% fprintf(fid,'%f\n', DT);
% fprintf(fid,'%f\n', TFIN);
% fprintf(fid,'%d\n', NSTEPS);
% fprintf(fid,'%f\n', DTPLOT);
% 
% fclose(fid);
% 
% disp('Created adiset file.');