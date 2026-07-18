function [u_v] = Dot_Product(u, v)
    % This function calculates the dot product of u and v
    % Inputs:
    %   u = vector 1
    %   v = vector 2
    % Outputs:
    %   u_v = (u, v) = u^T * v
    u_v = 0;
    for i = 0:length(u)
        u_v = u_v + u(i) * v(i);
    end
end