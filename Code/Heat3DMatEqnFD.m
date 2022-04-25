function CFS = Heat3DMatEqnFD(f,N,T,alph)
% Size of the grid; N should be even so that r=0 is not a grid point.
N = N+mod(N,2); P = N+1; 
r = linspace(-1,1,N+2)'; r = r(2:N+1); hr = 2/(N+1);
z = r;
th = pi*trigpts(P);
[rr,zz,tt] = ndgrid(r,z,th);

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

% Store values of initial function in (r, z, th) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
CFS = zeros(size(zz));
finit = @(r,z,th) f(r,z,th,0);
for j=1:N
    CFS(:,j,:) = finit(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

% Fourier transform in theta to decouple theta and r variables.
CFS(:,:,:,1) = V2C_cyl(CFS,'00t');
CFS(:,:,:,2:T)=zeros;

%1/10 time step; BDF2 30 times!


Hist = cell(3,1);

% Plug in exact solution for the first four time steps.

for j = 1:3
    for k=1:N
        Hist{j}(:,k,:) = f(rr(:,k,:),zz(:,k,:), tt(:,k,:),j*alph);
    end
    CFS(:,:,:,j+1)=V2C_cyl(Hist{j},'rzt');
end


%t=5:T; full time step
q = (12/25)*alph;
for d = 5:T
        RHS = (48/25)*CFS(:,:,:,d-1)-(36/25)*CFS(:,:,:,d-2)+(16/25)*CFS(:,:,:,d-3)-(3/25)*CFS(:,:,:,d-4);
        CFS(:,:,:,d) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

% Convert back to value space
CFS=C2V_cyl(CFS,'00t0');
end

%%
function Fapprox = Heat3DSolver(cfs,Q,N,C02,A,D2,coF)
P = N + 1;
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
