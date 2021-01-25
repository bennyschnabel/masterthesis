function [M, radii] = mil_tensor(fileName)
% MIL_TENSOR Calculation of the MIL tensor according to the approach of 
% Z. Tabor in "On the equivalence of two methods of determining fabric 
% tensor" from the year 2009.
%
%   [M] = MIL_TENSOR(fileName)
%   fileName ... File name

[~, ~, fExt] = fileparts(fileName);

switch lower(fExt)
  case '.csv'
    importData = table2array(readtable(fileName));
    MIL = importData(:,1);
    theta = importData(:,2);
    phi = importData(:,3);
  case '.dat'
    importData = importdata(fileName);
    MIL = importData(:,3);
    theta = importData(:,1);
    phi = importData(:,2);
  otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
    error('Unexpected file extension: %s', fExt);
end

x = zeros(length(MIL), 1);
y = x;
z = x;

for kk = 1 : 1 : length(x)
    [x(kk), y(kk), z(kk)] = sc2cc(MIL(kk), theta(kk), phi(kk));
end

[center, radii, evecs, v, chi2] = ellipsoid_fit( [ x y z ], '' );
fprintf( '\nAverage deviation of the fit: %.5f\n', sqrt( chi2 / size( x, 1 ) ) );

ae = v(1);
be = v(2);
ce = v(3);
de = v(4);
ee = v(5);
fe = v(6);

M = [[ae, de, ee]; [de, be, fe]; [ee, fe, ce]];

try chol(M);
    disp('MIL tensor M is symmetric positive definite.')
catch ME
    disp('MIL tensor M is not symmetric positive definite')
end
end