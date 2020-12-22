clc; clear all; close all;

% Implementation for calculating the MIL tensor in 3D
%
% Abbreviations used:
% I ... Image
% ROI ... Region Of Interest
% ROIBW ... Region Of Interest Black White (Binary)
%

%% User Input

% Image file name
imageFileName = 'Knochenprobe21.mat';
% Number of randomly generated angles (positiv integer)
numberOfDifferentAngles = 250;
% Distance between two created lines (positiv integer)
increment = 1;

%% Check if Matlab or GNU Octave

if isOctave() == 0
    % Matlab
    addpath('functions/');
    addpath('images/');
    addpath('matlab2tikz/');
elseif isOctave() == 1
    % GNU Octave
    addpath ("images:")
    addpath ("functions:")
else
    disp('Error')
    return
end

%% Import file


I = cell2mat(struct2cell(load(imageFileName)));
[r, c, p] = size(I);


%r = 10;
%c = 10;
%p = 10;

%% Test

x0 = [1; 0; 0];
y0 = [0; 1; 0];
z0 = [0; 0; 1];

P0 = [0; 0; 0];
P1 = [1; 0; 0];

theta = deg2rad(45);
alpha = deg2rad(-45);

P2 = rot3dx3(theta) * P1;
P2 = rot3dx2(alpha) * P2;

%P2 = [c; r; p];

[X,Y,Z] = bresenham_3d(P1, P2);

%% Plot
%{
plot3([P0(1) x0(1)], [P0(2) x0(2)], [P0(3) x0(3)], 'k')
hold on
plot3([P0(1) y0(1)], [P0(2) y0(2)], [P0(3) y0(3)], 'k')
hold on
plot3([P0(1) z0(1)], [P0(2) z0(2)], [P0(3) z0(3)], 'k')
hold on
plot3([P0(1) P1(1)], [P0(2) P1(2)], [P0(3) P1(3)], 'r-')
hold on
plot3([P0(1) P2(1)], [P0(2) P2(2)], [P0(3) P2(3)], 'r--')
xlabel('x_{1}')
ylabel('x_{2}')
zlabel('x_{3}')
axis equal
% [xmin xmax ymin ymax zmin zmax]
%axis([1 c 1 r 1 p])
%}
%% Test

n = 30;
m = 40;
o = 50;

%I = randi([0 1], n,m,o);

%[r, c, p] = size(I);

d = sqrt(r^2 + c^2 + p^2);

theta = deg2rad(90);
phi = deg2rad(0);
ra = 1;

P0 = [0; 0; 0];
[x, y, z] = sc2cc(ra, theta, phi);
P1 = [x; y; z];
n = P1 - P0;
n = round(1 / norm(n) * (n), 4);
[X,Y,Z] = bresenham_3d(P0, P1);
inc = 1;
[MIL] = calculate_mil_3d(n, r, c, p, inc, I);
disp(MIL)

figure(1)

xlabel('x_{1}')
ylabel('x_{2}')
zlabel('x_{3}')