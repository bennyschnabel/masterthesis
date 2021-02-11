function [BVTV] = volume_fraction(I)

% Bone volume: 0

BV = sum(I(:) == 0);
TV = sum(I(:) == 0) + sum(I(:) == 1);
BVTV = BV / TV;
end