% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

% This is the driver for the MacCormack solvers, it handles setting up the
% solver, running the solver, and saving the data.  This script also
% handles the variable speed portion of the problem.  This will work with 
% the 2-2 and 2-4 versions of MacCormack with the corresponding lines in
% the loop needing to be uncommented.

% n is the number of steps in the solution's domain
n = 201;

% The boundary conditions
xl = -1;
xr = 1;

% The parameters needed for the problem
sigma = 50; %sqrt(50.0);

% The step size for the solution grid
dx = (xr - xl) / (n - 1);

% I don't know what nzero does yet
nzero = (n - 1) / 2 + 1;

% This loop is going to populate the needed vectors for the structure
% of the solution
x = zeros(1, n);
c = zeros(1, n);
u = zeros(1, n);
uinit = zeros(1, n);

% j = n is the same as j = 1
for j = 1:n
    x(j) = xl + dx * (j - 1);
    % c(j) = 1.0;

    % for variable speed
    if (abs(x(j)) <= 0.25)
      c(j) = 1.0 + .9 * cos(2.0 * pi * x(j));
    else
      c(j) = 1.0;
    end

    arg = -(sigma * x(j))^2;
    u(j) = exp(arg);
    uinit(j) = exp(arg);
end

% These parameters are needed to find the solution
tbeg = 0;
tfin = 40.5;
cfl = .2;
cmax = 1.9;
dt = cfl * dx / cmax;
lambda = dt / dx;
dtplot = 10.0;
t = 0;
tplot = 0;
count = 0;

% This loop will solve the PDE using MacCormack's method which involves a
% step using a forwards-backwards step followed by a backwards-forward
% step.  In addition the solution will be saved periodically.
iplot = 0;
while (t < tfin)
    % REDACTED: solver code is available on request
end
hold off

plot(time,uzero)
for j = 1:n
    xplot0(j,1) = x(j);
    xplot0(j,2) = uinit(j);
    xplot1(j,1) = x(j);
    xplot1(j,2) = u10(j);
    xplot2(j,1) = x(j);
    xplot2(j,2) = u20(j);
    xplot3(j,1) = x(j);
    xplot3(j,2) = u30(j);
    xplot4(j,1) = x(j);
    xplot4(j,2) = u40(j);
end

for j = 1:count
    xplott(j,1) = time(j);
    xplott(j,2) = uzero(j);
end

save('plotx.dat', 'plotx', '-ascii');
save('uinit.dat', 'xplot0', '-ascii');
save('u10.dat', 'xplot1', '-ascii');
save('u20.dat', 'xplot2', '-ascii');
save('u30.dat', 'xplot3', '-ascii');
save('u40.dat', 'xplot4', '-ascii');
save('uzero.dat', 'xplott', '-ascii');