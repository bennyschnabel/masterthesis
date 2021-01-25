function show_ellipsoid(fileName)

[color1, color2, color3] = import_custom_colors();

% 
[~, ~, fExt] = fileparts(fileName);

switch lower(fExt)
  case '.csv'
    importData = table2array(readtable(fileName));
    MIL = importData(:,1);
    theta = importData(:,2);
    phi = importData(:,3);
  case '.dat'
    importData = importdata(fileName);
    MIL = importData(:,3);
    theta = importData(:,1);
    phi = importData(:,2);
  otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
    error('Unexpected file extension: %s', fExt);
end

x = zeros(length(MIL), 1);
y = x;
z = x;

for kk = 1 : 1 : length(x)
    [x(kk), y(kk), z(kk)] = sc2cc(MIL(kk), theta(kk), phi(kk));
end

[ ~, radii, evecs, v, chi2 ] = ellipsoid_fit( [ x y z ], '' );

%% Plot

figure()
for ii = 1 : 1 : length(x)
plot3([0 x(ii)], [0 y(ii)], [0 z(ii)], '.', 'Color', [color1(1), color2(1), color3(1)])
hold on
end

mind = min( [ x y z ] ) * 1.25;
maxd = max( [ x y z ] ) * 1.25;
nsteps = 100;
step = ( maxd - mind ) / nsteps;
[ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
hold off;
set( p, 'FaceColor', [color1(2), color2(2), color3(2)], 'EdgeColor', 'none' );
view(30,30)
grid on
alpha(0.3)   
axis vis3d equal;
camlight;
lighting phong;
xlabel('x_{1}')
ylabel('x_{2}')
zlabel('x_{3}')
end