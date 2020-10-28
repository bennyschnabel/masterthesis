function [M, ellipse, definite] = mean_intercept_length(fileName)

importData = table2array(readtable(fileName));

MIL = importData(:,1);
alpha = importData(:,3);

[alpha,order]=sort(alpha);
MIL=MIL(order);

% Convert to cartesian coordinates

x = zeros(length(importData), 1);
y = x;

for ii = 1 : length(importData)
    x(ii) = MIL(ii) * cos(alpha(ii));
    y(ii) = MIL(ii) * sin(alpha(ii));
end

ellipse = fit_ellipse(x, y);

a = ellipse(1);
b = ellipse(2);
c = ellipse(3);
d = ellipse(4);
e = ellipse(5);
f = ellipse(6);

%%

% Matrix of the quadratic equation
AQ = [[a, b/2, d/2]; [b/2, c, e/2]; [d/2, e/2, f]];

% Matrix of the quadratic form
A33 = [[a, b/2]; [b/2, c]];

% Calculate Eigenvectors and Eigenvalues
[~, eval] = eig(A33);

isEllipse = 0;

if det(AQ) ~=  0
    if det(A33) > 0
        disp('det(A33) > 0 -> Ellipse')
        isEllipse = 1;
    elseif det(A33) < 0
        disp('det(A33) < 0 -> Parabola')
    else
        disp('det(A33) < 0 -> hyperbola')
    end
else
    disp('The conic is degenerate')
end

if isEllipse == 1
    if eval(1,1) > 0 && eval(2,2) > 0
        disp('M is said to be positive definite')
        definite = 1;
    elseif eval(1,1) >= 0 && eval(2,2) >= 0
        disp('M is said to be positive semi-definite')
        definite = 2;
    elseif eval(1,1) < 0 && eval(2,2) < 0
        disp('M is said to be negative definite')
        definite = 3;
    elseif eval(1,1) <= 0 && eval(2,2) <= 0
        disp('M is said to be negative semi-definite')
        definite = 4;
    else
        disp('M is said to be indefinite')
        definite = 5;
    end
    M = A33;
else
    M = [[NaN, NaN]; [NaN, NaN]];
end
end