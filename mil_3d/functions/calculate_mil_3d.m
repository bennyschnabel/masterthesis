function [MIL] = calculate_mil_3d(n, r, c, p, inc, I)

% P = [r, c, p];
dispFigure = 'true';

Ctau = 0;
h = 0;

if round(n(1), 3) == 1 && round(n(2), 3) == 0 && round(n(3)) == 0
    disp('100')
    for pp = 1 : inc : p
        for pc = 1 : inc : c
            P1 = [1; pc; pp];
            fac = (r - P1(1)) / n(1);
            P2 = round(P1 + fac * n);
            
            [x,y,z] = bresenham_3d(P1, P2);
            
            h = h + norm(P2 - P1);
            
            switch dispFigure
                case 'true'
                    x_plot = x;
                    y_plot = zeros(size(x));
                    z_plot = zeros(size(x));
                    for ii = 1 : length(x)
                        y_plot(ii) = NaN;
                        z_plot(ii) = NaN;
                    end
                    
                    for ii = 1 : length(x) - 1
                        if I(x(ii), y(ii), z(ii)) == 0 && I(x(ii + 1), y(ii + 1), z(ii + 1)) == 1
                            Ctau = Ctau + 1;
                            x_plot(ii) = x(ii);
                            y_plot(ii) = y(ii);
                            z_plot(ii) = z(ii);
                        end
                    end
                    
                    figure(1)
                    plot3(x_plot, y_plot, z_plot, 'kx')
                    axis([0 r 0 c 0 p])
                    %pause(0.1)
                    hold on
                    
                case 'false'
                    for ii = 1 : length(x) - 1
                        if I(x(ii), y(ii), z(ii)) == 0 && I(x(ii + 1), y(ii + 1), z(ii + 1)) == 1
                            Ctau = Ctau + 1;
                        end
                    end
                otherwise
                    warning('Todo')
            end
        end
    end
elseif round(n(1), 3) == 0 && round(n(2), 3) == 1 && round(n(3)) == 0
    disp('010')
elseif round(n(1), 3) == 0 && round(n(2), 3) == 0 && round(n(3)) == 1
    disp('001')
elseif round(n(1), 3) == -1 && round(n(2), 3) == 0 && round(n(3)) == 0
    disp('-100')
end

MIL = h / Ctau;

end