function [imgmatrix, r, c] = img2matrix(type, filename, imagePath)

% Bresenham-Algorithm
% [imgmatrix, r, c] = img2matrix(type, filename, imagePath)
% Suppoerted application:
% type = 'matlab' or type = 'octave'
% Convert loaded image to matrix
% filePath only required for GNU Octave

if ~exist('imagePath','var')
    % third parameter does not exist, so default it to something
    imagePath = 0;
end

if type == 'matlab'
    imgdata = imread(filename);
elseif type == 'octave'
    image = file_in_path(imagePath, filename);
    imgdata = imread(image);
else
    disp('Software not supported')
end

[r, c, cc] = size(imgdata);

for ii = 1 : cc - 1
    if imgdata(:,:,ii) ~= imgdata(:,:,ii + 1)
        errorMessage = ['Color channel not equal'];
        error(errorMessage)
    end
end

imgmatrix = imgdata(:,:,1);

for ii = 1 : r
    for jj = 1 : c
        imgmatrix(ii, jj) = imgmatrix(ii, jj) / 255;
    end
end
end