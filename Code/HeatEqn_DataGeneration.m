%% Exact solution
w54 = 18.98013387517992112;  % 4th root of the J-Bessel function of order 5
w33 = 13.0152007216984344;   % 3rd root of the J-Bessel function of order 5
q1 = 3;                      % Oscillation in z.
q2 = 2;                      

% Exact solution
% u = @(r,z,th,t) 10*exp(-((pi*q1)^2 + w33^2)*t)*sin(q1*pi*z).*besselj(3,w33*r).*sin(3*th) + ...
%                 exp(-((pi*q2)^2 + w54^2)*t)*sin(q2*pi*z).*besselj(5,w54*r).*cos(5*th);
% u = @(r,z,th,t) 10*exp(-((pi*q1)^2 + w33^2)*t)*sin(q1*pi*z).*besselj(3,w33*r).*sin(3*th);

u = @(r,z,th,t) exactSolution(r,z,th,t);

% Use dt = 1e-3, because the solution decays very quickly.
dt = 1e-9;
T = 200;
T_FD = 6;

% Evaluate at a n-by-n-by-n grid at time t=T*dt
finit = @(r,z,th) u(r,z,th,0*dt);
f = @(r,z,th) u(r,z,th,(T-1)*dt);

Sp_disc = [7 11 15 19];
Sp_discFD = [61 81 101 121];

%% Test for harmonic.
n = 61;
Ops = Operators(n);
D0 = Ops.D0;
R0 = Ops.R0;

INIT = func2grid(finit,n,'rzt');
INIT_coeff = V2C_cyl(INIT,'rzt');
INIT_lap = zeros(n,n,n);
for j = -(n-1)/2:(n-1)/2
    INIT_lap(:,:,j+(n+1)/2) = R0^2*INIT_coeff(:,:,j+(n+1)/2)*(D0.')^2 + R0*D0*INIT_coeff(:,:,j+(n+1)/2) + ...
        R0^2*D0^2*INIT_coeff(:,:,j+(n+1)/2) - j^2*INIT_coeff(:,:,j+(n+1)/2);
end
INIT_lap = C2V_cyl(INIT_lap,'rzt');
CHECK = zeros(n,n,n);

for j = -(n-1)/2:(n-1)/2
    CHECK(:,:,j+(n+1)/2) = -((pi*q1)^2 + w33^2).*R0^2*INIT_coeff(:,:,j+(n+1)/2);
end
CHECK = C2V_cyl(CHECK,'rzt');
disp(norm(INIT_lap(:) - CHECK(:),inf)/norm(CHECK(:),inf))

%% New method
tic
TEST = Heat3DMatEqn_ExactSln(u,Sp_disc(1),T,dt);
TIME(1,1) = toc;
exact = func2grid(f,Sp_disc(1),'rzt');
ERROR(1,1) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(1)^3,1),inf)/norm(exact(:),inf);
disp('small')

% Higher Size.
tic
TEST = Heat3DMatEqn_ExactSln(u,Sp_disc(2),T,dt);
TIME(2,1) = toc;
exact = func2grid(f,Sp_disc(2),'rzt');
ERROR(2,1) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(2)^3,1),inf)/norm(exact(:),inf);
% ERROR(2,1) = norm(reshape(func2grid(f,81,'rzt') - TEST(:,:,:,T),81^3,1),inf);
disp('medium')

% Higher size.
tic
TEST = Heat3DMatEqn_ExactSln(u,Sp_disc(3),T,dt);
TIME(3,1) = toc;
exact = func2grid(f,Sp_disc(3),'rzt');
ERROR(3,1) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(3)^3,1),inf)/norm(exact(:),inf);
% ERROR(3,1) = norm(reshape(func2grid(f,101,'rzt') - TEST(:,:,:,T),101^3,1),inf);
disp('large')

% Higher size.
tic
TEST = Heat3DMatEqn_ExactSln(u,Sp_disc(4),T,dt);
TIME(4,1) = toc;
exact = func2grid(f,Sp_disc(4),'rzt');
ERROR(4,1) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(4)^3,1),inf)/norm(exact(:),inf);
% ERROR(4,1) = norm(reshape(func2grid(f,121,'rzt') - TEST(:,:,:,T),121^3,1),inf);
disp('very large')
%% Collocation
tic
TEST = Heat3DMatEqnColloc_ExactSln(u,Sp_disc(1),T,dt);
TIME(1,2) = toc;
exact = func2grid(f,Sp_disc(1),'rzt');
ERROR(1,2) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(1)^3,1),inf)/norm(exact(:),inf);
disp('small')

% Higher Size.
tic
TEST = Heat3DMatEqnColloc_ExactSln(u,Sp_disc(2),T,dt);
TIME(2,2) = toc;
exact = func2grid(f,Sp_disc(2),'rzt');
ERROR(2,2) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(2)^3,1),inf)/norm(exact(:),inf);
% ERROR(2,1) = norm(reshape(func2grid(f,81,'rzt') - TEST(:,:,:,T),81^3,1),inf);
disp('medium')

% Higher size.
tic
TEST = Heat3DMatEqnColloc_ExactSln(u,Sp_disc(3),T,dt);
TIME(3,2) = toc;
exact = func2grid(f,Sp_disc(3),'rzt');
ERROR(3,2) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(3)^3,1),inf)/norm(exact(:),inf);
% ERROR(3,1) = norm(reshape(func2grid(f,101,'rzt') - TEST(:,:,:,T),101^3,1),inf);
disp('large')

% Higher size.
tic
TEST = Heat3DMatEqnColloc_ExactSln(u,Sp_disc(4),T,dt);
TIME(4,2) = toc;
exact = func2grid(f,Sp_disc(4),'rzt');
ERROR(4,2) = norm(exact(:)-reshape(TEST(:,:,:,T),Sp_disc(4)^3,1),inf)/norm(exact(:),inf);
% ERROR(4,1) = norm(reshape(func2grid(f,121,'rzt') - TEST(:,:,:,T),121^3,1),inf);
disp('very large')
%% Finite Difference
n = Sp_discFD(1)-1;
r = linspace(-1,1,n+2)'; r = r(2:n+1);
[rr,zz,tt] = ndgrid(r,r,pi*trigpts(n+1));
FF = f(rr,zz,tt);

tic
TEST = Heat3DMatEqnFD_ExactSln(u,n,T_FD,dt);
TIME(1,3) = toc;
ERROR(1,3) = norm(reshape(FF - TEST(:,:,:,T_FD),n^2*(n+1),1), inf)/norm(FF(:),inf);
disp('small')

% Higher Size.
n = Sp_discFD(2)-1;
r = linspace(-1,1,n+2)'; r = r(2:n+1);
[rr,zz,tt] = ndgrid(r,r,pi*trigpts(n+1));
FF = f(rr,zz,tt);

tic
TEST = Heat3DMatEqnFD_ExactSln(u,n,T_FD,dt);
TIME(2,3) = toc;
ERROR(2,3) = norm(reshape(FF - TEST(:,:,:,T_FD),n^2*(n+1),1), inf)/norm(FF(:),inf);
disp('medium')

% Higher size.
n = Sp_discFD(3)-1;
r = linspace(-1,1,n+2)'; r = r(2:n+1);
[rr,zz,tt] = ndgrid(r,r,pi*trigpts(n+1));
FF = f(rr,zz,tt);

tic
TEST = Heat3DMatEqnFD_ExactSln(u,n,T_FD,dt);
TIME(3,3) = toc;
ERROR(3,3) = norm(reshape(FF - TEST(:,:,:,T_FD),n^2*(n+1),1), inf)/norm(FF(:),inf);
disp('large')

% Higher size.
n = Sp_discFD(4)-1;
r = linspace(-1,1,n+2)'; r = r(2:n+1);
[rr,zz,tt] = ndgrid(r,r,pi*trigpts(n+1));
FF = f(rr,zz,tt);

tic
TEST = Heat3DMatEqnFD_ExactSln(u,n,T_FD,dt);
TIME(4,3) = toc;
ERROR(4,3) = norm(reshape(FF - TEST(:,:,:,T_FD),n^2*(n+1),1), inf)/norm(FF(:),inf);
disp('very large')