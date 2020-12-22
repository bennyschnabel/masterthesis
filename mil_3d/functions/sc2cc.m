function [x, y, z] = sc2cc(r, theta, phi)

x = r * sin(theta) * cos(phi);
y = r * sin(theta) * sin(phi);
z = r * cos(theta);
end

% Spherical coordinate
% Cartesian coordinates