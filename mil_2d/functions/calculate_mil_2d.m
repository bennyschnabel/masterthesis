function [MIL] = calculate_mil_2d(n, r, c, xs, ys, inc, I_ROIBW)

n180 = rot2d(pi) * n;

Ctau = 0;
h = 0;

V1 = [round(xs(1)); round(ys(1))];
V2 = [round(xs(2)); round(ys(2))];
V3 = [round(xs(3)); round(ys(3))];
V4 = [round(xs(4)); round(ys(4))];

if n(1) > 0 && n(2) > 0
    [x1, y1, ~] = bresenham(V2, V1);
    [x2, y2, ~] = bresenham(V3, V4);

    for ii = 1 : inc : length(x1)
        P1 = [x1(ii); y1(ii)];
        P2 = [x2(ii); y2(ii)];
        
        if P1(1) < 1
            fac = (1 - P2(1)) / n180(1);
            P1 = round(P2 + fac * n180);
        end

        if P1(2) < 1
            fac = (1 - P2(2)) / n180(2);
            P1 = round(P2 + fac * n180);
        end
        
        if P2(1) > c
            fac = (c - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
        end

        if P2(2) > r
            fac = (r - P1(2)) / n(2);
            P2 = round(P1 + fac * n);
        end
        
        if P1(1) == 1 && P1(2) > r
            P1 = [1; r];
        end
        
        if P1(1) > c && P1(2) == 1
            P1 = [c; 1];
        end
        
        if P2(1) < 1 && P2(2) == r
            P2 = [1; r];
        end
        
        if P2(1) == c && P2(2) < 1
            P2 = [c; 1];
        end

        [x, y, ~] = bresenham(P1, P2);
        
        h = h + norm(P2 - P1);
        %h = h + length(x);
        
        % For plotting
        x_plot = x;
        y_plot = zeros(size(x));
        for jj = 1 : length(x)
            y_plot(jj) = NaN;
        end
        
        for jj = 1 : length(x) - 1
            if I_ROIBW(y(jj), x(jj)) == 0 && I_ROIBW(y(jj + 1), x(jj + 1)) == 1
                Ctau = Ctau + 1;
                %x_plot(jj) = x(jj);
                %y_plot(jj) = y(jj);
            end
        end
        
        %{
        figure(2)
        plot(x_plot, y_plot, 'kx')
        axis([0 c 0 r])
        hold on
        pause(0.1)
        %}
    end
    
elseif n(1) < 0 && n(2) > 0
    [x1, y1, ~] = bresenham(V1, V4);
    [x2, y2, ~] = bresenham(V2, V3);
    
    for ii = 1 : inc : length(x1)
        P1 = [x1(ii); y1(ii)];
        P2 = [x2(ii); y2(ii)];
        
        if P1(1) > c
            fac = (c - P2(1)) / n(1);
            P1 = round(P2 + fac * n);
        end
        
        if P1(2) < 1
            fac = (1 - P2(2)) / n(2);
            P1 = round(P2 + fac * n);
        end
        
        if P2(1) < 1
            fac = (1 - P1(1)) / n180(1);
            P2 = round(P1 + fac * n180);
        end
        
        if P2(2) > r
            fac = (r - P1(2)) / n180(2);
            P2 = round(P1 + fac * n180);
        end
        
        if P1(1) < 1 && P1(2) == 1
            P1 = [1; 1];
        end
        
        if P1(1) == c && P1(2) > r
            P1 = [c; r];
        end
        
        if P2(1) == 1 && P2(2) < 1
            P2 = [1; 1];
        end
        
        if P2(1) > c && P2(2) == r
            P2 = [c; r];
        end

        [x, y, ~] = bresenham(P1, P2);
        
        h = h + norm(P2 - P1);
        %h = h + length(x);
        
        for jj = 1 : length(x) - 1
            if I_ROIBW(y(jj), x(jj)) == 0 && I_ROIBW(y(jj + 1), x(jj + 1)) == 1
                Ctau = Ctau + 1;
            end
        end
    end
    
elseif round(n(1)) == 1 && n(2) == 0
    for py = 1 : inc : r
        P1 = [1; py];
        fac = (c - P1(1)) / n(1);
        P2 = round(P1 + fac * n);

        [x, y, ~] = bresenham(P1, P2);
        
        h = h + norm(P2 - P1);
        %h = h + length(x);
        
        for jj = 1 : length(x) - 1
            if I_ROIBW(y(jj), x(jj)) == 0 && I_ROIBW(y(jj + 1), x(jj + 1)) == 1
                Ctau = Ctau + 1;
            end
        end
    end
    
% n = [0 1]
elseif n(1) == 0 && round(n(2)) == 1
    for px = 1 : inc : c
        P1 = [px; 1];
        fac = (r - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        [x, y, ~] = bresenham(P1, P2);
        
        h = h + norm(P2 - P1);
        %h = h + length(x);
        
        for jj = 1 : length(x) - 1
            if I_ROIBW(y(jj), x(jj)) == 0 && I_ROIBW(y(jj + 1), x(jj + 1)) == 1
                Ctau = Ctau + 1;
            end
        end
    end
    
elseif round(n(1)) == -1 && n(2) == 0
    for py = 1 : inc : r
        P1 = [c; py];
        fac = (1 - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        [x, y, ~] = bresenham(P1, P2);

        h = h + norm(P2 - P1);
        %h = h + length(x);

        for jj = 1 : length(x) - 1
            if I_ROIBW(y(jj), x(jj)) == 0 && I_ROIBW(y(jj + 1), x(jj + 1)) == 1
                Ctau = Ctau + 1;
            end
        end
    end
end

MIL = h / Ctau;

end