clc; clear all; close all;

%% User Input

% Image file name
imageFileName = 'Knochenprobe2.int1.stream.tiff';
% DPI value of the image
dpi = 1814;
% Side length of the square/ cube created from the image [mm]
length = 5;
% Select start image from stack (eg. '2' for second image in stack)
startImage = 1;
% Select if 2D or 3D [2d, 3d]
plotType = '3d';
% Resolution of the export [dpi]
resolution = 300;

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

%% Create or load stack

% Find dot in file name
dotLocations = find(imageFileName == '.');
fileNameCropped = imageFileName(1:dotLocations(1) - 1);

fileName = [fileNameCropped, '.mat'];

% Check if data-file exists
if isfile(fileName)
    load(fileName);
else
    img2array(imageFileName, fileName);
end

%% Calculate lengths

pixelX = size(tiff_stack, 1);
pixelY = size(tiff_stack, 2);

lengthX = pixelX * 25.4 / dpi;
lengthY = pixelY * 25.4 / dpi;
lengthZ = lengthX / pixelX;

%% Convert to binary image

values = histcounts(tiff_stack, 256);
level = otsu(values);

level = 1 / 256 * level;
I = imbinarize(tiff_stack, level);

x = pixelX / 2;
va = round(pixelX/lengthX * length);
xl = x - round(va / 2);
xu = x + round(va / 2);

% Export image

switch plotType
    case {'2d', '2D'}
        % Export 2D image

        I_ROIBW = I(xl:xu,xl:xu,startImage);

        imageFileName2d = [fileNameCropped, '_2d_', num2str(length), ...
            'mm_', num2str(startImage), '.tex'];

        figure()
        imagesc(I_ROIBW)
        titleString = ['Threshold \tau: ', num2str(level * 256)];
        %titleString = 'Initial image';
        title(titleString)
        xlabel('x_{1}')
        ylabel('x_{2}')
        colorbar
        colormap('gray')
        
        exportgraphics(gcf, imageFileName2d, 'Resolution',resolution)
        %matlab2tikz(imageFileName2d);
        
    case {'3d', '3D'}
        % Export 3D image

        I_ROIBW = I(xl:xu,xl:xu,startImage:va);

        imageFileName3d = [fileNameCropped, '_3d_', num2str(length), ...
            'mm_', num2str(startImage), '.png'];

        figure()
        patch(isosurface(I_ROIBW,0.5), 'FaceColor', [17 17 17]/255, ...
            'EdgeColor', 'none')
        xlabel('x_{1}')
        ylabel('x_{2}')
        zlabel('x_{3}')
        view(30,30)
        axis vis3d
        alpha(0.3)

        exportgraphics(gcf, imageFileName3d, 'Resolution',resolution)
    otherwise
        warning('Unexpected plot type. No plot created.')
end