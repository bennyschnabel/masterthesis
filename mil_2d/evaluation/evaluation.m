clc; clear all; close all;

addpath('../')

T = readtable('export_mil.csv');

% exportData = [numberOfRuns, executionTime, ellipse(1), ellipse(2), ...
%    ellipse(3), ellipse(4), ellipse(5), ellipse(6)];
A = table2array(T);

numberOfRuns = [100];%, 250];

dexecutionTime = A(:,2);
da = A(:,3);
db = A(:,4);
dc = A(:,5);
dd = A(:,6);
de = A(:,7);
df = A(:,8);
ddefinite = A(:,9);

% Definition of colors
color1 = [83/255; 0/255; 228/255];
color2 = [103/255; 190/255; 175/255];
color3 = [108/255; 255/255; 52/255];

for jj = 1 : 1 : length(numberOfRuns)

[row, col] = find(A(:,1) == numberOfRuns(jj));

executionTime = zeros(length(row), 1);
a = zeros(length(row), 1);
b = zeros(length(row), 1);
c = zeros(length(row), 1);
d = zeros(length(row), 1);
e = zeros(length(row), 1);
f = zeros(length(row), 1);

for ii = 1 : 1 : length(row)
    executionTime(ii) = dexecutionTime(row(ii));
    a(ii) = da(row(ii));
    b(ii) = db(row(ii));
    c(ii) = dc(row(ii));
    d(ii) = dd(row(ii));
    e(ii) = de(row(ii));
    f(ii) = df(row(ii));
end

aSorted = sort(abs(a));
bSorted = sort(abs(b));
cSorted = sort(abs(c));

figure(1)
mu = arithmetic_mean(aSorted);
var = variance(aSorted, mu);
sigma = sqrt(var);
f = probability_density_function(aSorted, sigma, mu);
displayNameString = ['a_{', num2str(numberOfRuns(jj)), '}: ', num2str(length(row)), ', \mu: ', num2str(mu)];
plot(aSorted, f, '-', 'color', [color1(jj), color2(jj), color3(jj)], 'LineWidth', 2, 'DisplayName', displayNameString)
hold on
xlabel('Parameter a')
ylabel('\phi_{\mu, \sigma^{2}} (t)')
legend('Location', 'Best')
grid on
%matlab2tikz('evaluation_a.tex')
end