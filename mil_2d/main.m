clc; clear all; close all;

% Implementation for calculating the MIL tensor in 2D
%
% Abbreviations used:
% I ... Image
% ROI ... Region Of Interest
% ROIBW ... Region Of Interest Black White (Binary)
%

%% User Input

% Image file name
imageFileName = 'knochenprobe_1.png';
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

%%

fileName = [imageFileName(1:end-4), '_', num2str(numberOfDifferentAngles), '.csv'];

fileNameExport = ['export_', imageFileName(1:end-4), '_', ...
    num2str(numberOfDifferentAngles), '.csv'];

I = imread(imageFileName);
I_ROI = I(:,:,1);

values = histcounts(I, 256);
level = otsu(values);
level = 1 / 256 * level;
I_ROIBW = imbinarize(I_ROI, level);

%show_histogram(I_ROI, level)

%{
figure()
imagesc(I_ROIBW)
titleString = ['Threshold value \tau: ', num2str(level * 256)];
title(titleString)
xlabel('x_{1}')
ylabel('x_{2}')
colorbar
colormap('gray')
set(gca,'YDir','normal')
%}

[r, c] = size(I_ROIBW);

%% Loop to calculate the value MIL(theta) 

for kk = 1 : 1 : numberOfDifferentAngles
    
    P1 = [1; 0];
    % Creation of a random angle theta in the interval [0, pi]
    theta = 0 + (pi - 0) .* rand(1, 1);
    [R] = rot2d(theta);
    P2 = R * P1;
    % Creation of a unit vector in direction of theta
    n = round(1 / norm(P2) * (P2), 4);

    % Creation of the corner points for the framework
    [xs, ys] = generate_corner_points(n, r, c);
    
    [MIL] = calculate_mil_2d(n, r, c, xs, ys, increment, I_ROIBW);
    
    dispString = ['kk: ', num2str(kk), '/', num2str(numberOfDifferentAngles), ...
        ', tau = ', num2str(round(rad2deg(theta), 1)), ', MIL = ', num2str(MIL)];
    disp(dispString)
    
    exportData = [MIL, theta];
    dlmwrite(fileName, exportData, '-append');
end

%% Create ellipse

[beta1, beta2, phi] = ellipse_equation(fileName);
[M] = mil_tensor(beta1, beta2, phi);
[H] = fabric_tensor(M);

%{
if numberOfDifferentAngles == 50
    delete knochenprobe_1_50.csv
elseif numberOfDifferentAngles == 100
    delete knochenprobe_1_100.csv
elseif numberOfDifferentAngles == 250
    delete knochenprobe_1_250.csv
elseif numberOfDifferentAngles == 500
    delete knochenprobe_1_500.csv
end
%}

[v, e] = eig(M);
[v1, e1] = eig(H);
exportData = [v(1,1), v(2,1), v(1,2), v(2,2), e(1,1), e(2,2), ...
    v1(1,1), v1(2,1), v1(1,2), v1(2,2), e1(1,1), e1(2,2)];
dlmwrite(fileNameExport, exportData, '-append');

show_ellipse(fileName, M, beta1, beta2, phi)

%delete *.csv

%% Further investigation

alpha_w = sum(I_ROIBW(:) == 1);

alpha_b = sum(I_ROIBW(:) == 0);

[rho_app] = average_apparent_density(alpha_w, alpha_b);

[nu] = poissons_coefficient(alpha_b, alpha_w);