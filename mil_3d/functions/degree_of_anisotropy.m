function [DA] = degree_of_anisotropy(r)
% DEGREE_OF_ANISOTROPY todo
% https://imagej.net/BoneJ
%
%   [DA] = degree_of_anisotropy(r)
%   r ... radii of ellipsoid as vector (e.g. r = [r1; r2; r3])

a = min(r);
c = max(r);

DA = 1 - ((1 / c^2) / (1 / a^2));
%DA = 1 - (a / c);
end