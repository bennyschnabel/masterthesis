function img2array(imageFileName, fileName)
tiff_info = imfinfo(imageFileName);
tiff_stack = imread(imageFileName, 1);

% Concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(imageFileName, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

save(fileName, 'tiff_stack')
end