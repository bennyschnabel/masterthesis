function [beta1, beta2, phi] = ellipse_equation(fileName)

importData = table2array(readtable(fileName));

MIL = importData(:,1);
tau = importData(:,2);

[tau,order]=sort(tau);
MIL=MIL(order);

% Convert to cartesian coordinates
     
x = zeros(length(importData), 1);
y = x;
xN = zeros(length(x) * 2, 1);
yN = xN;

for ii = 1 : length(importData)
    x(ii) = MIL(ii) * cos(tau(ii));
    xN(ii) = x(ii);
    y(ii) = MIL(ii) * sin(tau(ii));
    yN(ii) = y(ii);
end

l = length(importData);
for ii = 1 : length(importData)
    xN(l + ii) = -x(ii);
    yN(l + ii) = -y(ii);
end

%%

XY = [xN, yN];
% Fit an ellipse using the Bookstein constraint
[z, beta1, beta2, phi] = fitellipse(XY, 'linear');
[z, beta1, beta2, phi] = fitellipse(XY);

x0 = z(1);
y0 = z(2);

[EL] = EllipseFitByTaubin(XY);
A = EL(1);
B = EL(2);
C = EL(3);
D = EL(4);
E = EL(5);
F = EL(6);
 
[beta1,beta2,x0,y0,phi] = f2(A,B,C,D,E,F);


[A,B,C,D,E,F] = f1(beta1,beta2,x0,y0,phi);

% Matrix of the quadratic form
A33 = [[A, B/2]; [B/2, C]];

AQ = [[A, B/2, D/2]; [B/2, C, E/2]; [D/2, E/2, F]];

%% Classification

if det(AQ) == 0
    disp('The conic is degenerate')
elseif det(AQ) ~= 0
    if det(A33) < 0
        disp('Q is a hyperbola')
    elseif det(A33) == 0
        disp('Q is a parabola')
    else
        if A == C && B == 0
            disp('Q is a circle')
        else
            if (A + C) * det(AQ) < 0
                disp('Q is an ellipse (real)')
            elseif (A + C) * det(AQ) > 0
                disp('Q is an ellipse (imaginary )')
            end
        end
    end
end

%% Cholesky-Zerlegung

try chol(A33);
    disp('Matrix A33 is symmetric positive definite.')
    positivDefinit = 1;
catch ME
    disp('Matrix A33 is not symmetric positive definite')
    positivDefinit = 0;
end

function [A,B,C,D,E,F] = f1(a,b,x0,y0,alpha)
A = a^2 * sin(alpha)^2 + b^2 * cos(alpha)^2;
B = 2 * (b^2 - a^2) * sin(alpha) * cos(alpha);
C = a^2 * cos(alpha)^2 + b^2 * sin(alpha)^2;
D = -2 * A * x0 - B * y0;
E = -B * x0 - 2 * C * y0;
F = A * x0^2 + B * x0 * y0 + C * y0^2 - a^2 * b^2;
end

function [a,b,x0,y0,theta] = f2(A,B,C,D,E,F)
a = (-sqrt(2 * (A*E^2 + C*D^2 - B*D*E + (B^2 - 4*A*C) * F) * ((A + C) + sqrt((A - C)^2 + B^2)))) / (B^2 - 4 * A * C);
b = (-sqrt(2 * (A*E^2 + C*D^2 - B*D*E + (B^2 - 4*A*C) * F) * ((A + C) - sqrt((A - C)^2 + B^2)))) / (B^2 - 4 * A * C);

x0 = (2 * C * D - B * E) /(B^2 - 4 * A * C);
y0 = (2 * A * E - B * D) /(B^2 - 4 * A * C);

if B ~= 0
    theta = atan(1 / B * (C - A - sqrt((A - C)^2 + B^2)));
elseif B == 0 && A < C
    theta = 0;
elseif B == 0 && A > C
    theta = pi/2;
end
end
end