function [alpha] = calculate_angle(n)

% Calculate angle
% [alpha] = calculate_angle()
% Return angle in rad

if n(1) > 0 && n(2) > 0
    alpha = acos(n(1) / norm(n)) * 180 / pi;
elseif n(1) < 0 && n(2) > 0
    alpha = acos(n(1) / norm(n)) * 180 / pi;
elseif n(1) < 0 && n(2) < 0
    alpha = 180 + (180 - (acos(n(1) / norm(n)) * 180 / pi));
elseif n(1) > 0 && n(2) < 0
    alpha = 270 + (90 - (acos(n(1) / norm(n)) * 180 / pi));
elseif round(n(1)) == 1.0 && n(2) == 0
    alpha = acos(n(1) / norm(n)) * 180 / pi;
elseif n(1) == 0 && round(n(2)) == 1.0
    alpha = acos(n(1) / norm(n)) * 180 / pi;
elseif round(n(1)) == -1.0 && n(2) == 0
    alpha = acos(n(1) / norm(n)) * 180 / pi;
elseif n(1) == 0 && round(n(2)) == -1.0
    alpha = 180 + acos(n(1) / norm(n)) * 180 / pi;
end

end