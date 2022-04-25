%%
% Size of the grid; N should be even so that r=0 is not a grid point.
N = 50; P = N+mod(N+1,2); T = 15; alph = .01;
r = linspace(-1,1,N+2)'; r = r(2:N+1); hr = 2/(N+1);
z = r;
th = pi*trigpts(P);
[zz,rr,tt] = meshgrid(z,r,th);

% Second order accurate differentiation matrices for d/dr and d^2/dr^2
D = (1/2/hr)*spdiags(ones(N,1)*[-1 1],[-1 1],N,N);
D2 = (1/hr/hr)*spdiags(ones(N,1)*[1 -2 1],-1:1,N,N);

% Multiplication matrices in value space by R and R^2
R = diag(r); R2 = R.^2;
I = speye(N);

C02 = I;
A=R2*D2+R*D;
D2 = D2;
coF=R2;

% Test Problem
f = @(r,z,t) r.^3.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

% Store values of initial function in (r, z, th) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
FF = zeros(size(zz));
for j=1:N
    FF(:,j,:) = f(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

% Fourier transform in theta to decouple theta and r variables.
CFS=zeros(N,N,P,T);
CFS(:,:,:,1) = V2C_cyl(FF,'00t');
Hist = zeros(N,N,P,30);

%1/10 time step; BDF2 30 times!


RHS = CFS(:,:,:,1);
q = alph/(10*T);
Hist(:,:,:,1) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);%Heat3DSolver(cfs,Q,N,C02,A,D2,coF)


for p = 1:29
        RHS = Hist(:,:,:,p);
        Hist(:,:,:,p+1) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

CFS(:,:,:,[2 3 4])=Hist(:,:,:,[10 20 30]);

%t=5:T; full time step
q = (12/25)*(alph/T);
for d = 5:T
        RHS = (48/25)*CFS(:,:,:,d-1)-(36/25)*CFS(:,:,:,d-2)+(16/25)*CFS(:,:,:,d-3)-(3/25)*CFS(:,:,:,d-4);
        CFS(:,:,:,d) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

% Convert back to value space
kapproxl=C2V_cyl(CFS,'00t0');
TimePlotEqui_cyl(kapproxl,.1, alph,.5);
%end

%%
function Fapprox = Heat3DSolver(cfs,Q,N,C02,A,D2,coF)
P = N + mod(N+1,2);
% Solves the screened Poisson equation, [I-Q*Del^2]u=CFS in coeff. space
% Homogeneous Dirichlet boundary conditions.

% Discretization size: N

Fapprox = zeros(N,N,P);
Ops = Operators(N,'DCT');
DCT = Ops.DCT;

bc = [ 1 zeros(1,N-1) ; 
       zeros(1,N-1) 1 ]; 

for k = -(P-1)/2 : (P-1)/2
    
    AA = coF - Q*(A - k^2*C02);
    BB = C02; 
    CC = -Q*coF; 
    DD = D2; 
    RHS = coF*cfs(:,:,k+(P+1)/2)*C02.';
    
    %RHS = RHS(1:N-2, 1:N-2);
    AA = AA(2:end-1,2:end-1);
    BB = BB(2:end-1,2:end-1);
    CC = CC(2:end-1,2:end-1);
    DD = DD(2:end-1,2:end-1);
    RHS = RHS(2:end-1,2:end-1);
    scale = norm(full(CC));
    CC = .2*CC/scale;
    DD = DD*scale*5;
    X22 = bartelsStewart(AA, BB, CC, DD, RHS );
    X = [ zeros(1,N) ; zeros(N-2,1) X22 zeros(N-2,1); zeros(1,N)]; 
    Fapprox(:, :, k+(P+1)/2) = X;

end

end
