function [MIL] = calculate_mil_3d(n, r, c, p, inc, I)
% CALCULATE_MIL_3D Calculation of the MIL values
%
%   [MIL] = CALCULATE_MIL_3D(n, r, c, p, inc, I)
%   MIL ... 
%   n ...
%   r ...
%   c ...
%   p ...
%   inc ...
%   I ...w

% P = [r, c, p];
dispFigure = 'false';

Ctau = 0;
h = 0;

if round(n(1), 3) == 1 && round(n(2), 3) == 0 && round(n(3), 3) == 0
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
                    
                    figure(99)
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
elseif round(n(1), 3) == 0 && round(n(2), 3) == 1 && round(n(3), 3) == 0
    disp('010')
    for pp = 1 : inc : p
        for pr = 1 : inc : r
            P1 = [pr; 1; pp];
            fac = (c - P1(2)) / n(2);
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
                    
                    figure(99)
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
elseif round(n(1), 3) == 0 && round(n(2), 3) == 0 && round(n(3), 3) == 1
    disp('001')
    for pc = 1 : inc : c
        for pr = 1 : inc : r
            P1 = [pr; pc; 1];
            fac = (p - P1(3)) / n(3);
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
                    
                    figure(99)
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
elseif round(n(1), 3) == -1 && round(n(2), 3) == 0 && round(n(3), 3) == 0
    disp('-100')
    for pp = 1 : inc : p
        for pc = 1 : inc : c
            P1 = [r; pc; pp];
            fac = (1 - P1(1)) / n(1);
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
                    
                    figure(99)
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
elseif round(n(1), 3) == 0 && round(n(2), 3) == -1 && round(n(3), 3) == 0
    disp('0-10')
    for pp = 1 : inc : p
        for pr = 1 : inc : r
            P1 = [pr; c; pp];
            fac = (1 - P1(2)) / n(2);
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
                    
                    figure(99)
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
elseif round(n(1), 3) == 0 && round(n(2), 3) == 0 && round(n(3), 3) == -1
    disp('00-1')
    for pc = 1 : inc : c
        for pr = 1 : inc : r
            P1 = [pr; pc; p];
            fac = (1 - P1(3)) / n(3);
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
                    
                    figure(99)
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
else
    % Room diagonal
    dR = norm([1;c;p] - [r;1;1]);
    % Radius of the sphere (Half room diagonal plus one)
    dR2 = dR / 2 + 1;
    % Origin of the sphere (center of the read data)
    PM = [round(r/2); round(c/2); round(p/2)];
    % Point on the sphere in direction n
    PS = round(PM + dR2 * n);
    
    % Create 4 straight lines orthogonal to the direction vector n,
    % spanning the plane Q and O, respectively.
    
    % TODO scale of the directional vector
    fac = dR2 / sqrt(2);
    % First directional vector (direction 1, 0)
    dV = [1; 0; (-n(1) * 1 - n(2) * 0) / (n(3))];
    % Normalize directional vector
    dV = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
    % First corner point of the plane Q
    Q1 = round(PS + fac * dV);
    % Second directional vector (direction 0, 1)
    dV = [0; 1; (-n(1) * 0 - n(2) * 1) / (n(3))];
    % Normalize directional vector
    dV = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
    % Second corner point of the plane Q
    Q2 = round(PS + fac * dV);
    % Third directional vector (direction 0, -1)
    dV = [0; -1; (-n(1) * 0 - n(2) * -1) / (n(3))];
    % Normalize directional vector
    dV = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
    % Third corner point of the plane Q
    Q3 = round(PS + fac * dV);
    % Fourth directional vector (direction -1, 0)
    dV = [-1; 0; (-n(1) * -1 - n(2) * 0) / (n(3))];
    % Normalize directional vector
    dV = (1 / abs(sqrt(dV(1)^2 + dV(2)^2 + dV(3)^2))) * dV;
    % Fourth corner point of the plane Q
    Q4 = round(PS + fac * dV);
    
    % TODO scale of the directional vector
    fac = dR2 * 2;
    % Counter vector to direction vector n
    n180 = -n;
    % First corner point of the plane O
    O1 = round(Q1 + fac * n180);
    % Second corner point of the plane O
    O2 = round(Q2 + fac * n180);
    % Third corner point of the plane O
    O3 = round(Q3 + fac * n180);
    % Fourth corner point of the plane O
    O4 = round(Q4 + fac * n180);
    
    % Generate meshgrid
    [x1Q12,x2Q12,x3Q12] = bresenham_3d(Q1, Q2);
    [x1Q34,x2Q34,x3Q34] = bresenham_3d(Q3, Q4);
    [x1O12,x2O12,x3O12] = bresenham_3d(O1, O2);
    [x1O34,x2O34,x3O34] = bresenham_3d(O3, O4);
    
    for kk = 1 : 1 : length(x1Q12)
        q12 = [x1Q12(kk); x2Q12(kk); x3Q12(kk)];
        q34 = [x1Q34(kk); x2Q34(kk); x3Q34(kk)];
        o12 = [x1O12(kk); x2O12(kk); x3O12(kk)];
        o34 = [x1O34(kk); x2O34(kk); x3O34(kk)];
        
        [x1q,x2q,x3q] = bresenham_3d(q12, q34);
        [x1o,x2o,x3o] = bresenham_3d(o12, o34);
        
        for ll = 1 : 1 : length(x1q)
            q = [x1q(ll); x2q(ll); x3q(ll)];
            o = [x1o(ll); x2o(ll); x3o(ll)];
            
            [x1, x2, x3] = bresenham_3d(q, o);
            counter = 0;
            indizes = [];
            for mm = 1 : inc : length(x1)
                if (x1(mm) >=  1 && x1(mm) <= r) && (x2(mm) >=  1 && x2(mm) <= c) ...
                    && (x3(mm) >=  1 && x3(mm) <= p)
                    counter = counter + 1;
                    indizes = [indizes, mm];
                end
            end
            
            if counter > 2
                SP = [x1(indizes(1)); x2(indizes(1)); x3(indizes(1))];
                EP = [x1(indizes(end)); x2(indizes(end)); x3(indizes(end))];
                
                h = h + norm(EP - SP);
                
                [x,y,z] = bresenham_3d(SP, EP);
                
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

                        figure(99)
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
    end
end

switch dispFigure
    case 'true'
        xlabel('x_{1}')
        ylabel('x_{2}')
        zlabel('x_{3}')
        grid on
        
        exportgraphics(gcf, 'mil_3d_bw_grid_4.png', 'Resolution',300)
    case 'false'
        disp('')
    otherwise
        warning('Todo')
end

MIL = h / Ctau;

end