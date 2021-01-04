clc; clear all; close all;

% Image file name
imageFileName = 'Knochenprobe2_1mm_1.mat';

I = cell2mat(struct2cell(load(imageFileName)));
[r, c, p] = size(I);


%outputFileName = 'Knochenprobe2_5mm_1.tif';
for K=1:p
    outputFileName = ['Knochenprobe2_5mm_1_', num2str(K,'%03d'),'.tif'];
    FILENAME = fullfile('bonej/', outputFileName);
    %imwrite(I(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
    imwrite(I(:, :, K), FILENAME, 'Compression','none');
end

