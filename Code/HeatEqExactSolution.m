%% Exact solution
w54 = 18.98013387517992112;  % 4th root of the J-Bessel function of order 5
w33 = 13.0152007216984344;   % 3rd root of the J-Bessel function of order 5
q1 = 3;                      % Oscillation in z.
q2 = 2;                      

% Exact solution
u = @(r,z,th,t) 10*exp(-((pi*q1)^2 + w33^2)*t)*sin(q1*pi*z).*besselj(3,w33*r).*sin(3*th) + ...
                exp(-((pi*q2)^2 + w54^2)*t)*sin(q2*pi*z).*besselj(5,w54*r).*cos(5*th);

% Use dt = 1e-3, because the solution decays very quickly.
dt = 1e-3;

% Evaluate at a n-by-n-by-n grid at time t=20*dt
n = 101;
[rr,zz,tt] = ndgrid(chebpts(n),chebpts(n),pi*trigpts(n));
uu = u(rr,zz,tt,dt*20);


%% Plot solution
VAL = uu;
N = size(VAL,1);
t = [pi*trigpts( N ); pi];
r = chebpts( N );
z = chebpts( N );
r = r(ceil(N/2):N);
[rr, zz, tt] = ndgrid(r,z,t);

VAL = VAL(ceil(N/2):N, :, :);
VAL = cat(3,VAL, VAL(:,:,1));

rslice = rr(floor(N/4)+1,1,1);
zslice = squeeze(zz(1,[floor(N/4)+1 floor(N/2)+1],1));
tslice = squeeze(tt(1,1,[1 floor(N/4)+1  floor(3*N/4)+1]+2)); 

hslicer = slice(zz,rr,tt,VAL,zslice,rslice,tslice);
hold on;
for j = 1:numel(hslicer)
        h = hslicer(j);
        [xs,ys,zs] = pol2cart(h.ZData,h.YData,h.XData);
        surf(xs,ys,zs,h.CData,'EdgeColor','none','FaceColor','Interp');
        hold on
end
delete(hslicer);
axis([-1 1 -1 1 -1 1]);
daspect([1 1 1])
hold off
% Add lighting and shine
camlight; lighting phong
% Change the colormap
colormap(jet)
colorbar;
% caxis([-scale scale]);
