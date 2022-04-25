% function p = verticalCylinderPatch(a,b,n,zb,zt)
function p = verticalCylinderPatch(xspan,yspan,zspan,n,clr)
thse = atan2(yspan,xspan);
% Handle the case where the patch crosses the negative x-axis.
if abs(diff(thse)) > pi
    thse = [thse(1) thse(1)+2*(pi-thse(1))];
end
ts = linspace(thse(1),thse(2),n)'; rs = 1+0*ts;
xt = rs.*cos(ts); yt = rs.*sin(ts);
p.x = [xt;flipud(xt)];
p.y = [yt;flipud(yt)];
p.z = [zspan(1)+0*xt;zspan(2)+0*yt];
p.th = [ts;flipud(ts)];
p.r = [ones(size(ts));ones(size(ts))];
p.clr = clr;

end
