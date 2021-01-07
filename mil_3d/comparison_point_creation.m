clc; clear all; close all;

%%

ru = 5;

theta = deg2rad(180*rand(1,1));
phi = deg2rad(360*rand(1,1));
ra = 1;

P0 = [0; 0; 0];
[x, y, z] = sc2cc(ra, theta, phi);
P1 = [x; y; z];
n = round(1 / norm(P1 - P0) * (P1 - P0), 4);

m = [0; 0; 0];

gs = m + ru * n;

%%

dV = [1; 0; (-n(1) * 1 - n(2) * 0) / (n(3))];
dV1 = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
Q1 = round(gs + (ru/2) * dV1);

dV = [0; 1; (-n(1) * 0 - n(2) * 1) / (n(3))];
dV2 = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
Q2 = round(gs + (ru/2) * dV2);

dV = [0; -1; (-n(1) * 0 - n(2) * -1) / (n(3))];
dV3 = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
Q3 = round(gs + (ru/2) * dV3);

dV = [-1; 0; (-n(1) * -1 - n(2) * 0) / (n(3))];
dV4 = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
Q4 = round(gs + (ru/2) * dV4);

[phi] = angle2vectors(n, dV1);
disp(phi)
[phi] = angle2vectors(n, dV2);
disp(phi)
[phi] = angle2vectors(n, dV3);
disp(phi)
[phi] = angle2vectors(n, dV4);
disp(phi)

%%

%{

plot3([m(1) gs(1)], [m(2) gs(2)], [m(3) gs(3)], '-k')
hold on
pause(0.5)
plot3([gs(1) Q1(1)], [gs(2) Q1(2)], [gs(3) Q1(3)], '-r')
hold on
pause(0.5)
plot3([gs(1) Q2(1)], [gs(2) Q2(2)], [gs(3) Q2(3)], '-g')
hold on
pause(0.5)
plot3([gs(1) Q3(1)], [gs(2) Q3(2)], [gs(3) Q3(3)], '-b')
hold on
pause(0.5)
plot3([gs(1) Q4(1)], [gs(2) Q4(2)], [gs(3) Q4(3)], '-c')
%}

% G : a + l * b
% gs = m + ru * n;
[M] = rot3axis(m, n, pi/2);
v = [Q1(1); Q1(2); Q1(3); 1];
Q2 = M * v;
Q2 = [Q2(1); Q2(2); Q2(3)];

[M] = rot(m, n, pi);
Q3 = M * v;
Q3 = [Q3(1); Q3(2); Q3(3)];

[M] = rot(m, n, (3/2) * pi);
Q4 = M * v;
Q4 = [Q4(1); Q4(2); Q4(3)];


plot3([m(1) gs(1)], [m(2) gs(2)], [m(3) gs(3)], '-k')
hold on
pause(0.5)
plot3([gs(1) Q1(1)], [gs(2) Q1(2)], [gs(3) Q1(3)], '-r')
hold on
pause(0.5)
plot3([gs(1) Q2(1)], [gs(2) Q2(2)], [gs(3) Q2(3)], '-g')
hold on
pause(0.5)
plot3([gs(1) Q3(1)], [gs(2) Q3(2)], [gs(3) Q3(3)], '-b')
hold on
pause(0.5)
plot3([gs(1) Q4(1)], [gs(2) Q4(2)], [gs(3) Q4(3)], '-c')



function [M] = rot(a, b, phi)

% G : a + l * b

d = sqrt(b(1)^2 + b(2)^2);

M = T(a) * Rzt(d, b) * Ryp(d, b) * Rz(phi) * Rypm(d, b) * Rztm(d, b) * T(-a);

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

function [phi] = angle2vectors(u, v)

phi = acos((dot(u, v)) / (norm(u) * norm(v)));
phi = rad2deg(phi);
end