function [R] = rot3axis(a, b, phi)
% ROT3AXIS Rotation around any axis
% Rotation R for the rotation of a point P around an arbitrarily oriented 
% axis G in space by an angle phi.
% [R] = rot3axis(a, b, phi)
% G : a + lambda * b
% <https://www.informatik.uni-leipzig.de/bsv/homepage/sites/default/files/CG_2.1-2.3_1.pdf>

d = sqrt(b(1)^2 + b(2)^2);

R = T(a) * Rzt(d, b) * Ryp(d, b) * Rz(phi) * Rypm(d, b) * Rztm(d, b) * T(-a);

function [T] = T(a)
T = [[1, 0, 0, a(1)]; [0, 1, 0, a(2)]; [0, 0, 1, a(3)]; [0, 0, 0, 1]];
end

function [R] = Rzt(d, b)
R = (1 / d) * [[b(1), -b(2), 0, 0]; [b(2), b(1), 0, 0]; [0, 0, d, 0]; [0, 0, 0, d]];
end

function [R] = Rztm(d, b)
R = (1 / d) * [[b(1), b(2), 0, 0]; [-b(2), b(1), 0, 0]; [0, 0, d, 0]; [0, 0, 0, d]];
end

function [R] = Ryp(d, b)
R = [[b(3), 0, d, 0]; [0, 1, 0, 0]; [-d, 0, b(3), 0]; [0, 0, 0, 1]];
end

function [R] = Rypm(d, b)
R = [[b(3), 0, -d, 0]; [0, 1, 0, 0]; [d, 0, b(3), 0]; [0, 0, 0, 1]];
end

function [R] = Rz(phi)
R = [[cos(phi), -sin(phi), 0, 0]; [sin(phi), cos(phi), 0, 0,]; [0, 0, 1, 0]; [0, 0, 0, 1]];
end
end