function TimePlotEqui_cyl(VAL,scale,alph,delay)
% This program takes four inputs, in this order: a 3 or 4 dimensional
% array of values, the desired axis scale, the time-step to show the user,
% and the desired delay between frames. It plots a scalar function on the
% cylinder using these parameters; if the array is 3D, only the first two
% arguments are needed, and no change in time will be involved.
%
% August 25, 2017. Grady Wright, David Darrow.

SIZE = size(VAL);
if ( size(SIZE,2) > 4 || size(SIZE,2) < 3)
    error('Wrong number of dimensions in array.')
end
N = SIZE(1);
P = SIZE(3);
if ( size(SIZE,2) == 3 )
    VAL = permute(VAL,[1 3 2]);
    VAL(:,:,:,2) = VAL;
    T = 2;
    alph = 0;
    delay = 0;
else
    VAL = permute(VAL,[1 3 2 4]);
    T = SIZE(4);
end


t = pi*trigpts( P );
r = linspace(-1,1,N);
z = r;
[tt, rr, zz] = meshgrid(t, r, z);

% Remove the doubled up data in the grid
tt2 = tt(floor(N/2)+1:end,:,:);
rr2 = rr(floor(N/2)+1:end,:,:);
zz2 = zz(floor(N/2)+1:end,:,:);

% Slices in the cylinder to plot
rslice = rr2(floor(N/4)+1,1,1);
zslice = squeeze(zz2(1,1,[floor(N/4)+1 floor(N/2)+1]));
tslice = tt2(1,[1 floor(P/4)+1  floor(3*P/4)+1]+2,1)'; 

% Add in the theta=pi to close the cylinder in the plots
size(tt2)
tt2 = [tt2 pi*ones(ceil(N/2),1,N)];
rr2 = [rr2 rr2(:,1,:)];
zz2 = [zz2 zz2(:,1,:)];

for d = 1:T
    ff2 = reshape(VAL(:,:,:,d),[N P N]);
    % Remove the doubled up data in the grid
    ff2 = ff2(floor(N/2)+1:end,:,:);
    ff2 = [ff2 repmat(ff2(:,1),[1 1 N])];
    hslicer = slice(tt2,rr2,zz2,ff2,tslice,rslice,zslice);
    hold on
    for j = 1:numel(hslicer)
        h = hslicer(j);
        [xs,ys,zs] = pol2cart(h.XData,h.YData,h.ZData);
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
    title(sprintf('time = %1.3f',d*alph))
    colorbar;
    caxis([-scale scale]);
    pause(delay);
end