function [X,Y,Z] = bresenham_3d(P1, P2, precision)
   if ~exist('precision','var') || isempty(precision) || round(precision) == 0
      precision = 0;
      P1 = round(P1);
      P2 = round(P2);
   else
      precision = round(precision);
      P1 = round(P1*(10^precision));
      P2 = round(P2*(10^precision));
   end
   d = max(abs(P2-P1)+1);
   X = zeros(1, d);
   Y = zeros(1, d);
   Z = zeros(1, d);
   x1 = P1(1);
   y1 = P1(2);
   z1 = P1(3);
   x2 = P2(1);
   y2 = P2(2);
   z2 = P2(3);
   dx = x2 - x1;
   dy = y2 - y1;
   dz = z2 - z1;
   ax = abs(dx)*2;
   ay = abs(dy)*2;
   az = abs(dz)*2;
   sx = sign(dx);
   sy = sign(dy);
   sz = sign(dz);
   x = x1;
   y = y1;
   z = z1;
   idx = 1;
   if(ax>=max(ay,az))			% x dominant
      yd = ay - ax/2;
      zd = az - ax/2;
      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;
         if(x == x2)		% end
            break;
         end
         if(yd >= 0)		% move along y
            y = y + sy;
            yd = yd - ax;
         end
         if(zd >= 0)		% move along z
            z = z + sz;
            zd = zd - ax;
         end
         x  = x  + sx;		% move along x
         yd = yd + ay;
         zd = zd + az;
      end
   elseif(ay>=max(ax,az))		% y dominant
      xd = ax - ay/2;
      zd = az - ay/2;
      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;
         if(y == y2)		% end
            break;
         end
         if(xd >= 0)		% move along x
            x = x + sx;
            xd = xd - ay;
         end
         if(zd >= 0)		% move along z
            z = z + sz;
            zd = zd - ay;
         end
         y  = y  + sy;		% move along y
         xd = xd + ax;
         zd = zd + az;
      end
   elseif(az>=max(ax,ay))		% z dominant
      xd = ax - az/2;
      yd = ay - az/2;
      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;
         if(z == z2)		% end
            break;
         end
         if(xd >= 0)		% move along x
            x = x + sx;
            xd = xd - az;
         end
         if(yd >= 0)		% move along y
            y = y + sy;
            yd = yd - az;
         end
         z  = z  + sz;		% move along z
         xd = xd + ax;
         yd = yd + ay;
      end
   end
   if precision ~= 0
      X = X/(10^precision);
      Y = Y/(10^precision);
      Z = Z/(10^precision);
   end
   return;					% bresenham_line3d
