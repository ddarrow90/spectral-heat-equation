% Set up and plot the top of the Cylinder in Cartesian coordinates.
N = 1001;
t = linspace(-pi,pi,N);
r = chebpts( 2*N-1 ); r = r(N:end);
z = linspace(0,2,4);
[tt, rr] = meshgrid(t, r);
xx = rr.*cos(tt);
yy = rr.*sin(tt);
% Function with zeros that match the brick pattern on the top.
ff = cos(1.5*pi*xx).*cos(1.5*pi*yy);
% Set the values of the function to different integers in the bricks so
% that we can assign colors as follows:
% 0 = red, 1 = green, 2 = blue, 3 = white, 4 = yellow.
clrmap = [[1 0 0];[0 1 0];[0 0 1];[1 1 1];[1 1 0]];
gg = ff;
gg( gg >= 0 ) = 1;
gg( gg < 0 ) = 0;
gg(yy >= 1/3 & gg < 1) = 2;
gg(yy <= -1/3 & gg < 1) = 3;
gg(xx >= 1/3 & yy <= -1/3 ) = 4;
gg(xx <= -1/3 & yy <= -1/3 ) = 4;
gg(xx <= -1/3 & yy >= 1/3 ) = 3;
gg(xx <= -1/3 & yy <= 1/3 & yy(:,:) >= -1/3) = 2;

% Plot the top of the cylinder
figure(1)
surf(xx,yy,z(end)+0*xx,gg);
hold on
axis equal
axis off
view(2);
shading interp;
colormap(clrmap)
m = 100;
% Plot lines outlining the bricks.
xl = ones(m,1)*[1/3 -1/3 1/3 -1/3];
yl = linspace(0+10*eps,1,m)'*[1 1 -1 -1]*sqrt(1-1/3^2);
xl = [xl yl];
yl = [yl xl(:,1:4)];
% plot3([1/3 1/3],[-1 1]*sqrt(1-1/3^2),[1 1],'k-','LineWidth',4)
% plot3(-[1/3 1/3],[-1 1]*sqrt(1-1/3^2),[1 1],'k-','LineWidth',4)
plot3(xl,yl,z(end)+0*xl,'k-','LineWidth',4)
% plot3([-1 1]*sqrt(1-1/3^2),[1/3 1/3],[1 1],'k-','LineWidth',4)
% plot3([-1 1]*sqrt(1-1/3^2),-[1/3 1/3],[1 1],'k-','LineWidth',4)
plot3(cos(t),sin(t),z(end)+0*t,'k-','LineWidth',4)
% export_fig -m2 -transparent -png 'CylindricalRubiksCartesian'
hold off

%% Plot the sides of the cylinder in Cartesian coordinates.
n = 100;
% The bricks touch the outer cylinder at (a,b), (b,a), etc.
a = (2*sqrt(2))/3; b = 1/3;
zspan = [[z(1) z(2)];[z(2) z(3)];[z(3) z(4)]];
% The starting and ending positions of the patches for the different bricks
% on the vertical sides.
xs = [[a a];[a b];[b -b];[-b -a];[-a -a];[-a -b];[-b b];[b a]];
ys = [[-b b];[b a];[a a];[a b];[b -b];[-b -a];[-a -a];[-a -b]]; 
ptch = [];

% The colors to assign.  There are 2*9 + 3*8 total colors and I have
% removed the ones already assigned.
clrs = {'r','r','r','r','r','r','r','r',...
        'g','g','g','g','g','g',...
        'b','b','b','b','b','b',...
        'y','y','y','y','y','y',...
        'w','w','w','w','w','w','w'};
rng(7122005);
% Select the colors at random.
hash = randperm(length(clrs));
cnt = 1;
for k = 1:8
    ptch{k,1} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(1,:),n,clrs{hash(cnt)});
    cnt = cnt + 1;
    ptch{k,2} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(2,:),n,clrs{hash(cnt)});
    cnt = cnt + 1;
    ptch{k,3} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(3,:),n,clrs{hash(cnt)});
    cnt = cnt + 1;
end
% Plot the patches.
figure(1), hold on
for k = 1:8
    fill3(ptch{k,1}.x,ptch{k,1}.y,ptch{k,1}.z,ptch{k,1}.clr,'EdgeColor','k','LineWidth',4)
    fill3(ptch{k,2}.x,ptch{k,2}.y,ptch{k,2}.z,ptch{k,2}.clr,'EdgeColor','k','LineWidth',4)
    fill3(ptch{k,3}.x,ptch{k,3}.y,ptch{k,3}.z,ptch{k,3}.clr,'EdgeColor','k','LineWidth',4)
end
daspect([1 1 1])
view(3)
hold off

%% Plot the top of the cylinder in cylindrical coordinates.
figure(4); hold on;
surf(tt,rr,z(end)+0*tt,gg);
view(2);
shading interp;
colormap([[1 0 0];[0 1 0];[0 0 1];[1 1 1];[1 1 0]])
hold on;
% Plot the outlines of the bricks in polar coordinates.
[tl,rl] = cart2pol(xl,yl);
plot3(tl,rl,z(end)+0*rl,'k-','LineWidth',4)
plot3(t,1+0*t,z(end)+0*t,'k-','LineWidth',4)
set(gca,'FontSize',14);
xlabel('$\theta$','Interpreter','Latex');
ylabel('$r$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi]);
set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{'0','1/2','1'});
set(gca,'ZTick',z);
set(gca,'ZTickLabel',{'0','1/3','2/3','1'});
daspect([1 1/2 1])
hold off
% export_fig -m2 -transparent -png 'CylindricalRubiksPolar'

%% Plot the sides of the cylinder in cylindrical coordiantes.
% The bricks touch the outer cylinder at (a,b), (b,a), etc.
a = (2*sqrt(2))/3; b = 1/3;
zspan = [[z(1) z(2)];[z(2) z(3)];[z(3) z(4)]];
% The starting and ending positions of the patches for the different bricks
% on the vertical sides.
xs = [[a a];[a b];[b -b];[-b -a];[-a -1];[-1 -a];[-a -b];[-b b];[b a]];
ys = [[-b b];[b a];[a a];[a b];[b+10*eps 0];[0-10*eps -b];[-b -a];[-a -a];[-a -b]]; 
ptch = [];
cnt = 1;
for k = 1:9
    if k == 6 % patch 5 and 6 cross the negative axis and must be dealt with separately.
        ptch{k,1} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(1,:),n,ptch{5,1}.clr);
        ptch{k,2} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(2,:),n,ptch{5,2}.clr);
        ptch{k,3} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(3,:),n,ptch{5,3}.clr);
    else
        ptch{k,1} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(1,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
        ptch{k,2} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(2,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
        ptch{k,3} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(3,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
    end
end

figure(4), hold on;
for k = 1:9
    if k == 5 || k == 6
        fill3(ptch{k,1}.th,ptch{k,1}.r,ptch{k,1}.z,ptch{k,1}.clr,'EdgeColor','none')
        plot3(ptch{k,1}.th(1:n),ptch{k,1}.r(1:n),ptch{k,1}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,1}.th(n+1:end),ptch{k,1}.r(n+1:end),ptch{k,1}.z(n+1:end),'k-','LineWidth',4)
        fill3(ptch{k,2}.th,ptch{k,2}.r,ptch{k,2}.z,ptch{k,2}.clr,'EdgeColor','none')
        plot3(ptch{k,2}.th(1:n),ptch{k,2}.r(1:n),ptch{k,2}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,2}.th(n+1:end),ptch{k,2}.r(n+1:end),ptch{k,2}.z(n+1:end),'k-','LineWidth',4)
        fill3(ptch{k,3}.th,ptch{k,3}.r,ptch{k,3}.z,ptch{k,3}.clr,'EdgeColor','none')
        plot3(ptch{k,3}.th(1:n),ptch{k,3}.r(1:n),ptch{k,3}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,3}.th(n+1:end),ptch{k,3}.r(n+1:end),ptch{k,3}.z(n+1:end),'k-','LineWidth',4)
    else
        fill3(ptch{k,1}.th,ptch{k,1}.r,ptch{k,1}.z,ptch{k,1}.clr,'EdgeColor','k','LineWidth',4)
        fill3(ptch{k,2}.th,ptch{k,2}.r,ptch{k,2}.z,ptch{k,2}.clr,'EdgeColor','k','LineWidth',4)
        fill3(ptch{k,3}.th,ptch{k,3}.r,ptch{k,3}.z,ptch{k,3}.clr,'EdgeColor','k','LineWidth',4)
    end
end
fill3(pi*[1 1 1 1],[1 0 0 1],[z(end) z(end) z(1) z(1)],'k','FaceColor',0.7*[1 1 1],'EdgeColor','none')
fill3(-pi*[1 1 1 1],[1 0 0 1],[z(end) z(end) z(1) z(1)],'k','FaceColor',0.7*[1 1 1],'EdgeColor','none')
plot3([pi pi],[1 0],[z(end-1) z(end-1)],'k-','LineWidth',4)
plot3([pi pi],[1 0],[z(end-2) z(end-2)],'k-','LineWidth',4)
plot3([pi pi],[1 0],[z(end-3) z(end-3)],'k-','LineWidth',4)
plot3([pi pi],[1/3 1/3],[z(end) z(end-1)],'k-')
plot3([pi pi],[1/3 1/3],[z(end-1) z(end-2)],'k-')
plot3([pi pi],[1/3 1/3],[z(end-2) z(end-3)],'k-')
plot3(-[pi pi],[1 0],[z(end-1) z(end-1)],'k-','LineWidth',4)
plot3(-[pi pi],[1 0],[z(end-2) z(end-2)],'k-','LineWidth',4)
plot3(-[pi pi],[1 0],[z(end-3) z(end-3)],'k-','LineWidth',4)
plot3(-[pi pi],[1/3 1/3],[z(end) z(end-1)],'k-')
plot3(-[pi pi],[1/3 1/3],[z(end-1) z(end-2)],'k-')
plot3(-[pi pi],[1/3 1/3],[z(end-2) z(end-3)],'k-')
daspect([1 1/2 1])
axis tight
view(140,30)
hold off

%% Plot the top of the doubled up cylinder 
figure(3);
GG = gg;
% Double up the function on the top.
GG = [flipud(GG(2:N,(N+3)/2:end)) flipud(GG(2:N,1:(N+1)/2)); GG];
TT = [flipud(tt(2:N,:,1)); tt(:,:,1)];
RR = [-flipud(rr(2:N,:,1)); rr(:,:,1)];
surf(TT,RR,z(end)+0*TT,GG)
colormap(clrmap)
shading interp
hold on;
rlf = -rl; tlf = -tl;
plot3(tl,rl,z(end)+0*rl,'k-','LineWidth',4)
plot3(t,1+0*t,z(end)+0*t,'k-','LineWidth',4)
plot3(tlf,rlf,z(end)+0*rlf,'k-','LineWidth',4)
plot3(t,-1+0*t,z(end)+0*t,'k-','LineWidth',4)
set(gca,'FontSize',14);
xlabel('$\theta$','Interpreter','Latex');
ylabel('$r$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi]);
set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
set(gca,'YTickLabel',{'-1','-1/2','0','1/2','1'});
set(gca,'ZTick',z);
set(gca,'ZTickLabel',{'0','1/3','2/3','1'});
axis tight
daspect([1 1/2 1])
hold off
% export_fig -m2 -transparent -png 'CylindricalRubiksPolarDoubled'

%% Plot the vertical sides of the doubled up cylinder.
a = (2*sqrt(2))/3; b = 1/3;
xs = [[a a];[a b];[b -b];[-b -a];[-a -1];[-1 -a];[-a -b];[-b b];[b a]];
ys = [[-b b];[b a];[a a];[a b];[b+10*eps 0];[0-10*eps -b];[-b -a];[-a -a];[-a -b]]; 
ptch = [];
cnt = 1;
for k = 1:9
    if k == 6
        ptch{k,1} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(1,:),n,ptch{5,1}.clr);
        ptch{k,2} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(2,:),n,ptch{5,2}.clr);
        ptch{k,3} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(3,:),n,ptch{5,3}.clr);
    else
        ptch{k,1} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(1,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
        ptch{k,2} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(2,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
        ptch{k,3} = verticalCylinderPatch(xs(k,:),ys(k,:),zspan(3,:),n,clrs{hash(cnt)});
        cnt = cnt + 1;
    end
end

figure(3), hold on;
for k = 1:9
    if k == 5 || k == 6
        fill3(ptch{k,1}.th,ptch{k,1}.r,ptch{k,1}.z,ptch{k,1}.clr,'EdgeColor','none')
        plot3(ptch{k,1}.th(1:n),ptch{k,1}.r(1:n),ptch{k,1}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,1}.th(n+1:end),ptch{k,1}.r(n+1:end),ptch{k,1}.z(n+1:end),'k-','LineWidth',4)
        fill3(ptch{k,2}.th,ptch{k,2}.r,ptch{k,2}.z,ptch{k,2}.clr,'EdgeColor','none')
        plot3(ptch{k,2}.th(1:n),ptch{k,2}.r(1:n),ptch{k,2}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,2}.th(n+1:end),ptch{k,2}.r(n+1:end),ptch{k,2}.z(n+1:end),'k-','LineWidth',4)
        fill3(ptch{k,3}.th,ptch{k,3}.r,ptch{k,3}.z,ptch{k,3}.clr,'EdgeColor','none')
        plot3(ptch{k,3}.th(1:n),ptch{k,3}.r(1:n),ptch{k,3}.z(1:n),'k-','LineWidth',4)
        plot3(ptch{k,3}.th(n+1:end),ptch{k,3}.r(n+1:end),ptch{k,3}.z(n+1:end),'k-','LineWidth',4)
    else
        fill3(ptch{k,1}.th,ptch{k,1}.r,ptch{k,1}.z,ptch{k,1}.clr,'EdgeColor','k','LineWidth',4)
        fill3(ptch{k,2}.th,ptch{k,2}.r,ptch{k,2}.z,ptch{k,2}.clr,'EdgeColor','k','LineWidth',4)
        fill3(ptch{k,3}.th,ptch{k,3}.r,ptch{k,3}.z,ptch{k,3}.clr,'EdgeColor','k','LineWidth',4)
    end
end
gam = atan2(b,a);
fill3(pi*[1 1 1 1],[1 -1 -1 1],[z(end) z(end) z(1) z(1)],'k','FaceColor',0.7*[1 1 1],'EdgeColor','none')
fill3(-pi*[1 1 1 1],[1 -1 -1 1],[z(end) z(end) z(1) z(1)],'k','FaceColor',0.7*[1 1 1],'EdgeColor','none')
plot3([pi pi],[1 -1],[z(end-1) z(end-1)],'k-','LineWidth',4)
plot3([pi pi],[1 -1],[z(end-2) z(end-2)],'k-','LineWidth',4)
plot3([pi pi],[1 -1],[z(end-3) z(end-3)],'k-','LineWidth',4)
plot3([pi pi],[1/3 1/3],[z(end) z(end-1)],'k-')
plot3([pi pi],[1/3 1/3],[z(end-1) z(end-2)],'k-')
plot3([pi pi],[1/3 1/3],[z(end-2) z(end-3)],'k-')
plot3([pi pi],-[1/3 1/3],[z(end) z(end-1)],'k-')
plot3([pi pi],-[1/3 1/3],[z(end-1) z(end-2)],'k-')
plot3([pi pi],-[1/3 1/3],[z(end-2) z(end-3)],'k-')
plot3(-[pi pi],[1 -1],[z(end-1) z(end-1)],'k-','LineWidth',4)
plot3(-[pi pi],[1 -1],[z(end-2) z(end-2)],'k-','LineWidth',4)
plot3(-[pi pi],[1 -1],[z(end-3) z(end-3)],'k-','LineWidth',4)
plot3(-[pi pi],[1/3 1/3],[z(end) z(end-1)],'k-')
plot3(-[pi pi],[1/3 1/3],[z(end-1) z(end-2)],'k-')
plot3(-[pi pi],[1/3 1/3],[z(end-2) z(end-3)],'k-')
daspect([1 1/2 1])
% axis([-pi-0.05 pi+0.05 -1.05 1.05 -0.05 2.05])
view(140,30)
hold off