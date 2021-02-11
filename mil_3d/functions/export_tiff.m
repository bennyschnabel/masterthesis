clc; clear all; close all;

%% User Input

% Image file name
imageFileName = 'Knochenprobe2_1.2mm_1.mat';
% Export as single tiff files or tiff stream [single, stream]
exportType = 'stream';

%% Calc

fileName = [imageFileName(1:end-4), '_'];

I = cell2mat(struct2cell(load(imageFileName)));
[r, c, p] = size(I);

switch exportType
    case 'single'
        for K=1:p
            outputFileName = [fileName, num2str(K,'%03d'),'.tif'];
            FILENAME = fullfile('bonej/', outputFileName);
            imwrite(I(:, :, K), FILENAME, 'Compression','none');
        end
    case 'stream'
        for K=1:p
            outputFileName = [fileName, 'full.tif'];
            FILENAME = fullfile('bonej/', outputFileName);
            imwrite(I(:, :, K), FILENAME, 'WriteMode', 'append',  'Compression','none');
        end
    otherwise
        disp('other value')
end