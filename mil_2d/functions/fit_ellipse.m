function [a] = fit_ellipse(x, y)

% 
% NUMERICALLY STABLE DIRECT LEAST SQUARES FITTING OF ELLIPSES
% Radim Hal´?r? & Jan Flusser

% quadratic part of the design matrix
D1 = [x .^2, x .* y, y .^2];

% linear part of the design matrix
D2 = [x, y, ones(size(x))];

% quadratic part of the scatter matrix
S1 = D1' * D1;

% combined part of the scatter matrix
S2 = D1' * D2;

% linear part of the scatter matrix
S3 = D2' * D2;

% for getting a2 from a1
T = - inv(S3) * S2';

% reduced scatter matrix
M = S1 + S2 * T;

%premultiply by inv(C1)
M = [M(3, :) ./ 2; - M(2, :); M(1, :) ./ 2];

% solve eigensystem
[evec, eval] = eig(M);

% evaluate a’Ca
cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :) .^ 2;

% eigenvector for min. pos. eigenvalue
a1 = evec(:, find(cond > 0));

% ellipse coefficients
a = [a1; T * a1];

end