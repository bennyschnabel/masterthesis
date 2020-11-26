clc; clear all; close all;

%
% ROI ... Region Of Interest
% ODF ... Orientation-Dependent Feature
%

application = 'M';

if application == 'M'
    % Matlab
    addpath('functions/');
    addpath('images/');
    addpath('matlab2tikz/');
elseif application == 'O'
    % GNU Octave
    addpath ("images:")
    addpath ("functions:")
else
    disp('Error')
    return
end

numberOfDifferentAngles = 250;

imageFileName = 'image1_0.png';

fileNameDefinit = ['definit_', imageFileName(1:end-4), '_', num2str(numberOfDifferentAngles), '.csv'];

fileName = [imageFileName(1:end-4), '_', num2str(numberOfDifferentAngles), '.csv'];

I = imread(imageFileName);
I_ROI = I(:,:,1);

%show_histogram(I_ROI)
values = histcounts(I, 256);
level = otsu(values);

level = 1 / 256 * level;
I_ROIBW = imbinarize(I_ROI, level);

[r, c] = size(I_ROIBW);

for kk = 1 : 1 : numberOfDifferentAngles
    
    P1 = [1; 0];
    tau = 0 + (pi - 0) .* rand(1, 1);
    [R] = rot2d(tau);
    P2 = R * P1;
    n = round(1 / norm(P2) * (P2), 4);

    [xs, ys] = generate_corner_points(n, r, c);
    
    [MIL] = calculate_mil_2d(n, r, c, xs, ys, I_ROIBW);
    
    dispString = ['kk: ', num2str(kk), '/', num2str(numberOfDifferentAngles), ...
        ', tau = ', num2str(round(rad2deg(tau), 1)), ', MIL = ', num2str(MIL)];
    disp(dispString)
    
    exportData = [MIL, tau];
    dlmwrite(fileName, exportData, '-append');
end

%% Create ellipse

[AQ, positivDefinit, lambda1, lambda2, x0, y0, a, b, alpha] = ellipse_equation(fileName);

exportData = [positivDefinit, lambda1, lambda2, AQ(1,1), AQ(1,2), AQ(2,2), ...
    AQ(3,1), AQ(3,2), AQ(3,3)];
dlmwrite(fileNameDefinit, exportData, '-append');

show_ellipse(fileName, AQ, a, b, x0, y0, alpha)

%% Further investigation

alpha_w = sum(I_ROIBW(:) == 1);

alpha_b = sum(I_ROIBW(:) == 0);

[rho_app] = average_apparent_density(alpha_w, alpha_b);

[nu] = poissons_coefficient(alpha_b, alpha_w);