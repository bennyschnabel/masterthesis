function [x_out, y_out, m] = bresenham(P1, P2)

% Bresenham-Algorithm
% [x, y] = bresenham(P1, P2)
% with P = [x, y] representing a Pixel
% <a
% href="https://de.wikipedia.org/wiki/Bresenham-Algorithmus">Wikipedia</a>
% <a href="https://www.k-achilles.de/algorithmen/bresenham-gerade.pdf">Achilles</a>

% Input P1 and P2 hat to be an Array with 2 elements
validateattributes(P1,{'numeric'}, {'numel',2})
validateattributes(P2,{'numeric'}, {'numel',2})

if mod(P1(1), 1) ~= 0
    errorMessage = ['Input P1(1) = ', num2str(P1(1)), ', but has to be an integer!'];
    error(errorMessage)
elseif mod(P1(2), 1) ~= 0
    errorMessage = ['Input P1(2) = ', num2str(P1(2)), ', but has to be an integer!'];
    error(errorMessage)
elseif mod(P2(1), 1) ~= 0
    errorMessage = ['Input P2(1) = ', num2str(P2(1)), ', but has to be an integer!'];
    error(errorMessage)
elseif mod(P2(2), 1) ~= 0
    errorMessage = ['Input P2(2) = ', num2str(P2(2)), ', but has to be an integer!'];
    error(errorMessage)
end

dx = P2(1) - P1(1);
dy = P2(2) - P1(2);

m = dy / dx;

xstep = 1;
ystep = 1;

x = P1(1);
y = P1(2);

x_out = zeros(dx + 1, 1);
y_out = zeros(dx + 1, 1);
ii = 1;

if dx < 0
    dx = -dx;
    xstep = -1;
end

if dy < 0
    dy = -dy;
    ystep = -1;
end

a = 2 * dx;
b = 2 * dy;

if dy <= dx
    f = -dx;
    while x ~= P2(1)
        x_out(ii) = x;
        y_out(ii) = y;
        ii = ii + 1;
        f = f + b;
        if f > 0
            y = y + ystep;
            f = f - a;
        end
        x = x + xstep;
    end
else
    f = -dy;
    while y ~= P2(2)
        x_out(ii) = x;
        y_out(ii) = y;
        ii = ii + 1;
        f = f + a;
        if f > 0
            x = x + xstep;
            f = f - b;
        end
        y = y + ystep;
    end
end
x_out(ii) = x;
y_out(ii) = y;
end