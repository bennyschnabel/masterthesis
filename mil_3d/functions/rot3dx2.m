function [R] = rot3dx2(theta)
R = [[cos(theta) 0 sin(theta)]; [0 1 0];  [-sin(theta) 0 cos(theta)]];
end