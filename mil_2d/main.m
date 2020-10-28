clc; clear all; close all;
%% Implementation of the MIL tensor calculation in 2D
% ...

%% Input parameters

% type = 'matlab' or type = 'octave'
type = 'matlab';
% Add folder with example images
addpath('images/');
% String of the file name
importFileName = 'image3.png';
% Number of randomly generated direction vectors of length one
numberOfDirections = 10;
% Number of simulation runs
numberOfSimulations = 1;

% Create subfolder for results
if ~exist('results', 'dir')
    mkdir('results')
end
addpath('results/');

%% Main loop
for kk = 1 : 1 : numberOfSimulations
    disp('-------------------------------------')
    stringNumberOfSimulations = ['Current number of simulation runs: ', num2str(kk)];
    disp(stringNumberOfSimulations)

    %% Set timer
    start = tic;

    %% Load image as matrix

    [imgmatrix, r, c] = img2matrix(type, importFileName);

    %% Calculate MIL scalar, alpha, h, Cv

    date = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    exportDataFileName = [date, '_export_data', '.csv'];
    exportMILFileName = 'export_mil.csv';

    for ii = 1 : 1 : numberOfDirections

    % Generate a random direction vector of length one

    P1 = [randi(c); randi(r)];
    P2 = [randi(c); randi(r)];

    if P1(1) == P2(1) && P1(2) == P2(2)
        P1 = [randi(c); randi(r)];
        P2 = [randi(c); randi(r)];
    end

    n = 1 / norm(P2 - P1) * (P2 - P1);

    % Calculate angle alpha

    [alphaDeg] = calculate_angle(n);
    alphaRad = deg2rad(alphaDeg);

    % Calculate MIL scalar

    h = 0;
    Cv = 0;
    [h, Cv, ~] = get_direction(n, c, r, imgmatrix, h, Cv);

    MIL = h/Cv;

    % Export data to csv

    exportData = [MIL, alphaDeg, alphaRad, Cv, h, n(1), n(2)];
    dlmwrite(exportDataFileName, exportData, '-append');
    end

    %% Calculate MIL Tensor M and Ellipse from MIL scalar and alpha

    [M, ellipse, definite] = mean_intercept_length(exportDataFileName);

    executionTime = toc(start);

    exportData = [numberOfDirections, executionTime, ellipse(1), ellipse(2), ...
        ellipse(3), ellipse(4), ellipse(5), ellipse(6), definite];
    dlmwrite(exportMILFileName, exportData, '-append');
    movefile(exportDataFileName, 'results');

end