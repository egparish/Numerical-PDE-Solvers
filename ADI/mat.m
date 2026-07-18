% NOTE: Parts of the file are withheld but available on request.
% Structure, documentation, and descriptions are intact.  the file will
% not run as posted.

function [AX, BX, CX, AY, BY, CY] = mat(NX, NY, DT, D1, D2)
    % Using the input data this function calculates the A, B, and C bands
    % for the x and y implicit steps of ADI.  The values in these diagonal
    % bands come from the stencil used in 2D ADI.  In the x direction the
    % banded structure has the following coefficients
    %       [ ... 0 -d1dt/2h^2  1 + 2(d1dt/2h^2)  -d1dt/2h^2 0 ... ]
    % In the y direction the banded structure has the following
    % coefficients
    %       [ ... 0 -d2dt/2h^2 1 + 2(d2dt/2h^2) -d2dt/2h^2 0 ... ]

    % Inputs:
    %   NX = number of points in the x direction
    %   NY = number of points in the y direction
    %   D1 = first diffusivities for the x direction
    %   D2 = second diffusivities for the y direction
    %   DT = timestep
    % Outputs:
    %   AX = lower diagonal for the x implicit step
    %   BX = diagonal for the x implicit step
    %   CX = upper diagonal for the x implicit step
    %   AY = lower diagonal for the y implicit step
    %   BY = diagonal for the y implicit step
    %   CY = upper diagonal for the y implicit step

    % REDACTED: the initialization code is available on request

end