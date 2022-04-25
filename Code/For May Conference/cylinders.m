%% Cheb, not doubled
n = 15; 
r = chebpts(n,[0,1]); 
theta = linspace(0,2*pi,n);
[rr,tt] = meshgrid(r,theta);
z = chebpts(n);
[xx,yy] = pol2cart(tt,rr);
for j = 1:numel(z)
    plot3(xx,yy,z(j)*ones(size(xx)),'.k','markersize',10)
    hold on 
end
rr = ones(1000,1)';
tt = linspace(0,2*pi,1000);
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,1*ones(size(yy)),'k-', 'linewidth',2)
plot3(xx,yy,-1*ones(size(yy)),'k-', 'linewidth',2)
rr = ones(1000,1)';
tt = zeros(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
tt = pi*ones(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
axis equal 
axis off
view(0,20)
hold off

%% Cheb, doubled

n = 15; 
r = chebpts(n,[-1,1]); 
theta = linspace(0,2*pi,1);
[rr,tt] = meshgrid(r,theta);
z = linspace(-1,1,n);
[xx,yy] = pol2cart(tt,rr);
for j = 1:numel(z)
    plot3(xx,yy,z(j)*ones(size(xx)),'.b','markersize',10)
    hold on 
end
r = 0; 
theta = linspace(0,2*pi,1);
[rr,tt] = meshgrid(r,theta);
z = linspace(-1,1,n);
[xx,yy] = pol2cart(tt,rr);
for j = 1:numel(z)
    plot3(xx,yy,z(j)*ones(size(xx)),'.r','markersize',10)
    hold on 
end
rr = ones(1000,1)';
tt = linspace(0,2*pi,1000);
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,1*ones(size(yy)),'k-', 'linewidth',2)
plot3(xx,yy,-1*ones(size(yy)),'k-', 'linewidth',2)
rr = ones(1000,1)';
tt = zeros(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
tt = pi*ones(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
axis equal 
axis on
view(10,10)
hold off

%% Equispaced, not doubled

n = 23; 
r = linspace(0,1,n); 
theta = linspace(0,2*pi,n);
[rr,tt] = meshgrid(r,theta);
z = linspace(-1,1,n);
[xx,yy] = pol2cart(tt,rr);
for j = 1:numel(z)
    plot3(xx,yy,z(j)*ones(size(xx)),'.k','markersize',10)
    hold on 
end
rr = ones(1000,1)';
tt = linspace(0,2*pi,1000);
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,1*ones(size(yy)),'k-', 'linewidth',2)
plot3(xx,yy,-1*ones(size(yy)),'k-', 'linewidth',2)
rr = ones(1000,1)';
tt = zeros(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
tt = pi*ones(1000,1)';
[xx,yy] = pol2cart(tt,rr);
plot3(xx,yy,chebpts(1000)','k-', 'linewidth',2)
axis equal 
axis off
view(0,20)
hold off

%% Cross Section, not doubled

n = 20; 
r = chebpts(n,[0,1]); 
theta = linspace(0,2*pi,n);
[rr,tt] = meshgrid(r,theta);
[xx,yy] = pol2cart(tt,rr);
plot(xx,yy,'.k','markersize',10)
hold on 
rr = ones(1000,1)';
tt = linspace(0,2*pi,1000);
[xx,yy] = pol2cart(tt,rr);
plot(xx,yy,'k-', 'linewidth',2)
axis equal 
axis off
hold off

%% Cross Section, doubled

n = 20; 
r = chebpts(n,[-1,1]); 
theta = linspace(0,2*pi,n);
[rr,tt] = meshgrid(r,theta);
z = chebpts(n);
[xx,yy] = pol2cart(tt,rr);
plot(xx,yy,'.k','markersize',10)
hold on 
rr = ones(1000,1)';
tt = linspace(0,2*pi,1000);
[xx,yy] = pol2cart(tt,rr);
plot(xx,yy,'k-', 'linewidth',2)
axis equal 
axis off
hold off