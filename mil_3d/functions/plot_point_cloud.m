clc; clear all; close all;

%% User unput

fileName = 'Knochenprobe2_1mm_1_500_21_01_15_13_35_28.csv';

%%

importData = table2array(readtable(fileName));

[color1, color2, color3] = import_custom_colors();

MIL = importData(:,1);
theta = importData(:,2);
phi = importData(:,3);

x = zeros(length(MIL), 1);
y = x;
z = x;

for kk = 1 : 1 : length(x)
    [x(kk), y(kk), z(kk)] = sc2cc(MIL(kk), theta(kk), phi(kk));
end

figure()
plot3(x, y, z, '.', 'color', [color1(1), color2(1), color3(1)])
view(30,30)
grid on
axis vis3d equal;
camlight;
lighting phong;
xlabel('x_{1}')
ylabel('x_{2}')
zlabel('x_{3}')