clc; clear all; close all;

% Implementation for calculating the MIL tensor in 3D
%
% Abbreviations used:
% I ... Image
% ROI ... Region Of Interest
% ROIBW ... Region Of Interest Black White (Binary)
%

%% User Input

% Image file name ()
imageFileName = 'Knochenprobe2_1mm_1.mat';
% Number of randomly generated angles (positiv integer, minimum 9, )
numberOfDifferentAngles = 2;
% Distance between two created lines (positiv integer)
increment = 5;

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

%% File names creation

fileName = [imageFileName(1:end-4), '_', num2str(numberOfDifferentAngles), '.csv'];

%% Loop to calculate the value MIL(theta) 

for kk = 1 : 1 : numberOfDifferentAngles
    % Spherical coordinates, generate direction vector
    theta = deg2rad(180*rand(1,1));
    phi = deg2rad(360*rand(1,1));
    ra = 1;

    P0 = [0; 0; 0];
    [x, y, z] = sc2cc(ra, theta, phi);
    P1 = [x; y; z];
    n = round(1 / norm(P1 - P0) * (P1 - P0), 4);
    
    [MIL] = calculate_mil_3d(n, r, c, p, increment, I);
    
    dispString = ['kk: ', num2str(kk), '/', num2str(numberOfDifferentAngles), ...
        ', theta = ', num2str(round(rad2deg(theta), 1)), ...
        ', phi = ', num2str(round(rad2deg(phi), 1)), ...
        ', MIL = ', num2str(round(MIL, 1))];
    disp(dispString)
    
    exportData = [MIL, theta, phi];
    dlmwrite(fileName, exportData, '-append');
end

%% Calculate ellipsoid

[radii] = ellipse_equation(fileName);

%% Calculate MIL tensor M and fabric tensor H

%[M] = mil_tensor(beta1, beta2, phi);
%[H] = fabric_tensor(M);

%% Degree of anisotropy

[DA] = degree_of_anisotropy(radii);