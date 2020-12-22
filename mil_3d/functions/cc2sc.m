function [r, theta, phi] = cc2sc(x, y, z)

r = sqrt(x^2 + y^2 + z^2);
phi = atan(y / x);
theta = acos(z / r);
end

% Spherical coordinate
% Cartesian coordinates