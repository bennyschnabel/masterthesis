function [h, Cv, m] = get_direction(n, c, r, imgmatrix, h, Cv)

% n = [+ +]
if n(1) > 0 && n(2) > 0
    % X-Direction
    for px = 1 : 1 : c
        P1 = [px; 1];
        fac = (c - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        if P2(2) > r
            fac = (r - P1(2)) / n(2);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end

    % Y-Direction
    for py = 1 : 1 : r
        P1 = [1; py];
        fac = (r - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        if P2(1) > c
            fac = (c - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
% n = [- +]
elseif n(1) < 0 && n(2) > 0
    % X-Direction
    for px = c : -1 : 1
        P1 = [px; 1];
        fac = (1 - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        if P2(2) > r
            fac = (r - P1(2)) / n(2);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
    % Y-Direction
    for py = 1 : 1 : r
        P1 = [c; py];
        fac = (r - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        if P2(1) < 1
            fac = (1 - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
 
% n = [- -]
elseif n(1) < 0 && n(2) < 0
    % X-Direction
    for px = c : -1 : 1
        P1 = [px; r];
        fac = (1 - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        if P2(2) < 1
            fac = (1 - P1(2)) / n(2);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
    % Y-Direction
    for py = r : -1 : 1
        P1 = [c; py];
        fac = (1 - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        if P2(1) < 1
            fac = (1 - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
% n = [+ -]
elseif n(1) > 0 && n(2) < 0
    % X-Direction
    for px = 1 : 1 : c
        P1 = [px; r];
        fac = (c - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        if P2(2) < 1
            fac = (1 - P1(2)) / n(2);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
    % Y-Direction
    for py = r : -1 : 1
        P1 = [1; py];
        fac = (1 - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        if P2(1) > c
            fac = (c - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
        end
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
% n = [1 0]
elseif round(n(1)) == 1 && n(2) == 0
    for py = 1 : 1 : r
        P1 = [1; py];
        fac = (c - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
% n = [0 1]
elseif n(1) == 0 && round(n(2)) == 1
    for px = 1 : 1 : c
        P1 = [px; 1];
        fac = (r - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
    
% n = [-1 0]
elseif round(n(1)) == -1 && n(2) == 0
    for py = 1 : 1 : r
        P1 = [c; py];
        fac = (1 - P1(1)) / n(1);
        P2 = round(P1 + fac * n);
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end

% n = [0 -1]
elseif n(1) == 0 && round(n(2)) == -1
    for px = 1 : 1 : c
        P1 = [px; r];
        fac = (1 - P1(2)) / n(2);
        P2 = round(P1 + fac * n);
        
        h = h + norm(P2 - P1);
        
        [x, y, m] = bresenham(P1, P2);
        
        x_plot = x;
        y_plot = zeros(size(x));
        for ii = 1 : length(x)
            y_plot(ii) = NaN;
        end
        
        for ii = 1 : length(x) - 1
            if imgmatrix(y(ii), x(ii)) == 0 && imgmatrix(y(ii + 1), x(ii + 1)) == 1
                Cv = Cv + 1;
                x_plot(ii) = x(ii);
                y_plot(ii) = y(ii);
            end
        end
        
        %plot(x_plot, y_plot, 'kx')
        %hold on
    end
end
end