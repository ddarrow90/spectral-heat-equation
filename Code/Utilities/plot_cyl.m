function plot_cyl(VAL,scale,ord)
if (nargin < 3)
    ord = 'rzt';
end
VAL = C2V_cyl(VAL,ord);
N = size(VAL,1);
t = [pi*trigpts( N+1 ); pi];
r = chebpts( N );
z = chebpts( N );
r = r(ceil((N+1)/2):N);
r = [0; r]; % Append zero to avoid R crossing center line.
[rr, zz, tt] = ndgrid(r,z,t);

% Define zero value as average of neighbors.
VAL(ceil((N+1)/2)-1, :, :) = .5*VAL(ceil((N+1)/2)-1, :, :)+.5*VAL(ceil((N+1)/2), :, :);
VAL = VAL(ceil((N+1)/2)-1:N, :, :);
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
if nargin > 1
    caxis([-scale scale]);
end
end