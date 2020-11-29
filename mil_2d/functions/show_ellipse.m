function show_ellipse(fileName, M, beta1, beta2, phi)

[color1, color2, color3] = import_custom_colors();

importData = table2array(readtable(fileName));

MIL = importData(:,1);
theta = importData(:,2);

[theta,order] = sort(theta);
MIL = MIL(order);

% Convert to cartesian coordinates
     
x = zeros(length(importData), 1);
y = x;
xN = zeros(length(x) * 2, 1);
yN = xN;

for ii = 1 : length(importData)
    x(ii) = MIL(ii) * cos(theta(ii));
    xN(ii) = x(ii);
    y(ii) = MIL(ii) * sin(theta(ii));
    yN(ii) = y(ii);
end

l = length(importData);
for ii = 1 : length(importData)
    xN(l + ii) = -x(ii);
    yN(l + ii) = -y(ii);
end

[evec, ~] = eig(M);

v1 = evec(:,1) * (max(x) / 100 * 50);
v2 = evec(:,2) * (max(x) / 100 * 50);

figure()
dataPointsString = ['Datapoint: ', num2str(length(x))];
plot(xN, yN, '-', 'Color', [color1(1) color2(1) color3(1)], ...
    'DisplayName',dataPointsString)
hold on


t = linspace(0, 2*pi, 200);
x1 = beta1 * cos(phi) * cos(t) - beta1 * sin(phi) * sin(t);
y1 = beta1 * sin(phi) * cos(t) + beta2 * cos(phi) * sin(t);

displayNameString = ['\beta_{1}: ', num2str(round(beta1, 1)), ', \beta_{2}: ', ...
    num2str(round(beta2, 1)), ', \phi: ', num2str(round(rad2deg(phi), 1)), '°'];
plot(x1, y1, '-', 'Linewidth', 2, 'Color',  ...
    [color1(2) color2(2) color3(2)], 'DisplayName', displayNameString)
%hold on
%plot([0 v1(1)], [0 v1(2)], '-k', 'Linewidth', 1, 'DisplayName', 'EV_{1}')
%hold on
%plot([0 v2(1)], [0 v2(2)], '-', 'Color',  ...
%    [color1(3) color2(3) color3(3)], 'Linewidth', 1, 'DisplayName', 'EV_{2}')
xlabel('x_{1}')
ylabel('x_{2}')
legend()

axis equal
grid on

end