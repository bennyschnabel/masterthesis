function [PN] = rotPoint(v, P, P1, theta)
    % https://robotics.stackexchange.com/questions/12782/how-rotate-a-point-around-an-arbitrary-line-in-3d
    a = v(1);
    b = v(2);
    c = v(3);
    d = sqrt(b^2 + c^2);
    T = [[1, 0, 0, -P1(1)]; [0, 1, 0, -P1(2)]; [0, 0, 0, -P1(3)]; [0, 0, 0, 1]];
    Rx = [[1, 0, 0, 0]; [0, c/d, -b/d, 0]; [0, b/d, c/d, 0]; [0, 0, 0, 1]];
    Ry = [[d, 0, -a, 0]; [0, 1, 0, 0]; [a, 0, d, 0]; [0, 0, 0, 1]];
    Rz = [[cos(theta), -sin(theta), 0, 0]; [sin(theta), cos(theta), 0, 0]; [0, 0, 1, 0]; [0, 0, 0, 1]];
    
    P0 = [P(1); P(2); P(3); 1];
    
    PN = inv(T) * Rx.^-1 * Ry.^-1 * Rz * Ry * Rx * T * P0;
end