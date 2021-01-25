function [H] = fabric_tensor(M)
% FABRIC_TENSOR Calculation of the fabric tensor according to the approach of 
% Cowin, Stephen C in "Wolff's law of trabecular architecture at remodeling 
% equilibrium" from the year 1986.
%
%   [H] = FABRIC_TENSOR(M)
%   M ... MIL tensor

[v,D] = eig(inv(M));

D = sqrt(D);

H = v * D * v';

%H = M^(-1/2);

try chol(H);
    disp('Fabric tensor H is symmetric positive definite.')
catch ME
    disp('Fabric tensor H is not symmetric positive definite')
end
end