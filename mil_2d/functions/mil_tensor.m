function [M] = mil_tensor(a, b, phi)
% MIL_TENSOR Calculation of the MIL tensor according to the approach of 
% Z. Tabor in "On the equivalence of two methods of determining fabric 
% tensor" from the year 2009.
%
%   [M] = MIL_TENSOR(a, b, phi)
%   a ... semi-major axis
%   b ... semi-minor axis
%   phi ... rotation angle (the angle from the positive horizontal axis 
%           to the ellipse's major axis)

sin2 = 1/2 * (1 - cos(2 * phi));
cos2 = 1/2 * (1 + cos(2 * phi));

M11 = a^2 * cos2 + b^2 * sin2;
M12 = (a^2 - b^2) * cos(phi) * sin(phi);
M22 = a^2 * sin2 + a^2 * cos2;

M = [[M11, M12]; [M12, M22]];
end