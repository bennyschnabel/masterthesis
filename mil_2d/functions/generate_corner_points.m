function [xs, ys] = generate_corner_points(n, r, c)

% PD  ------  PC
%    |      |
%    |      |
%    |      |
% PA  ------  PB

t = 1;

PA = [t; t];
PB = [c; t];
PC = [c; r];
PD = [t; r];


xs = zeros(4, 1);
ys = zeros(4, 1);

% Check direction of unit vector

if n(1) < 0
    n = rot2d(-pi/2) * n;
end

% A - B

n180 = rot2d(pi) * n;
a1 = -n180(2);
b1 = n180(1);
c1 = PB(1) * a1 + PB(2) * b1;

n90 = rot2d(pi/2) * n;
a2 = -n90(2);
b2 = n90(1);
c2 = PA(1) * a2 + PA(2) * b2;

xs(1) = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
ys(1) = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

% A - D

a1 = -n180(2);
b1 = n180(1);
c1 = PD(1) * a1 + PD(2) * b1;

a2 = -n90(2);
b2 = n90(1);
c2 = PA(1) * a2 + PA(2) * b2;

xs(2) = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
ys(2) = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

% C - D

a1 = -n(2);
b1 = n(1);
c1 = PD(1) * a1 + PD(2) * b1;

a2 = -n90(2);
b2 = n90(1);
c2 = PC(1) * a2 + PC(2) * b2;

xs(3) = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
ys(3) = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

% B - C

a1 = -n(2);
b1 = n(1);
c1 = PB(1) * a1 + PB(2) * b1;

a2 = -n90(2);
b2 = n90(1);
c2 = PC(1) * a2 + PC(2) * b2;

xs(4) = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
ys(4) = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

end