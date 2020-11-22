function show_datapoints(fileName)

[color1, color2, color3] = import_custom_colors();

importData = table2array(readtable(fileName));

MIL = importData(:,1);
tau = importData(:,2);

[tau,order] = sort(tau);
MIL = MIL(order);

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

figure()
plot(x, y, 'x', 'Color', [color1(1) color2(1) color3(1)])
hold on 
plot(xN, yN, '-', 'Color', [color1(2) color2(2) color3(2)])
xlabel('x_{1}')
ylabel('x_{2}')
titleString = ['\Sigma Datenpunkte: ', num2str(length(x))];
title(titleString)
axis equal
grid on
matlab2tikz('sum_datapoints.tex')
end