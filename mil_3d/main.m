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
imageFileName = 'Knochenprobe2.mat';
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

%% 

fileName = [imageFileName(1:end-4), '_', num2str(numberOfDifferentAngles), '.csv'];

%% Loop to calculate the value MIL(theta) 

% Test values

n = 30;
m = 40;
o = 50;
%I = randi([0 1], n,m,o);
%[r, c, p] = size(I);

d = sqrt(r^2 + c^2 + p^2);
%{
theta = deg2rad(90);
phi = deg2rad(0);
ra = 1;

P0 = [0; 0; 0];
[x, y, z] = sc2cc(ra, theta, phi);
P1 = [x; y; z];
n = P1 - P0;
n = round(1 / norm(P1 - P0) * (P1 - P0), 4);
[X,Y,Z] = bresenham_3d(P0, P1);
inc = 1;
[MIL] = calculate_mil_3d(n, r, c, p, inc, I);
disp(MIL)
%}

%% Loop to calculate the value MIL(theta) 

angles = [0; 90; 180];
%angles = 180;
for kk = 1 : 1 : size(angles)
    
    theta = deg2rad(90);
    phi = deg2rad(angles(kk));
    ra = 1;

    P0 = [0; 0; 0];
    [x, y, z] = sc2cc(ra, theta, phi);
    P1 = [x; y; z];
    n = round(1 / norm(P1 - P0) * (P1 - P0), 4);
    
    [MIL] = calculate_mil_3d(n, r, c, p, increment, I);
    
    dispString = ['kk: ', num2str(kk), '/', num2str(size(angles)), ...
        ', theta = ', num2str(round(rad2deg(theta), 1)), ...
        ', phi = ', num2str(round(rad2deg(phi), 1)), ...
        ', MIL = ', num2str(MIL)];
    disp(dispString)
    
    exportData = [MIL, theta, phi];
    dlmwrite(fileName, exportData, '-append');
end

