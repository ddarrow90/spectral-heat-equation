% This program is a fast Stokes' Equation solver for the cylinder with 
% radius 1 and height 2. It takes 
% 
% and use a Fourier-Chebyshev-Chebyshev expansion. 
% 
% Author: David Darrow, February 2017.


% Cylindrical harmonic
%f = @(t,r,z) feval(diskfun.harmonic(4,3),t,r,'polar').*sin(2*pi*z);

% Discretization size, time discretization size: 
N = 9;
T = 500; % Number of time steps
alph = .01; % Time step size

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

% Construct grid (t,r,z) = Fourier-Chebyshev-Chebyshev:
t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );
[zz, rr, tt] = meshgrid(z, r, t);

ops = Operators(N,'C02 D2 A coF Dt iDCT iDFT');
C02 = ops.C02;
Dt = ops.Dt;
D2 = ops.D2;
A = ops.A;
coF = ops.coF;
iDCT = ops.iDCT;
iDFT = ops.iDFT;

H = HarmonicToRows(N,1e-14);


%%
% Initial (time 0) Values: 
tor = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
pol = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

% Store values of initial functions in (p, r, z) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
TOR = zeros(size(zz));
for j=1:N
    TOR(:,j,:) = tor(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

% Convert from value space to coefficient space
TOR = V2C_cyl(TOR,'rzt');
TOR(:,:,:,2:T)=zeros;

POL = zeros(size(zz));
for j=1:N
    POL(:,j,:) = pol(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

%TimePlot_cyl(POL(:,:,:,1),.05)

% Convert from value space to coefficient space
POL = V2C_cyl(POL,'rzt');
POL(:,:,:,2:T)=zeros;
% Project initial conditions onto nullspace of boundary matrix.
symm = .5*kron(spdiags((1-(-1).^(1:N)).',0,N,N),kron(speye(N),(speye(N)-flip(speye(N)))*iDCT));
symm = symm + .5*kron(spdiags((1+(-1).^(1:N)).',0,N,N),kron(speye(N),(speye(N)+flip(speye(N)))*iDCT));
symm = [symm; kron([ones(N-1,1) -speye(N-1)]*iDFT,kron(speye(N),[zeros(1,(N-1)/2) 1 zeros(1,(N-1)/2)]*iDCT))];
% symm: at even Fourier modes, F(r=1) = F(r=-1). At odd Fourier modes,
% F(r=1) = -F(r=-1). At r=0, F is constant in theta.
norm(symm*reshape(POL(:,:,:,1),N^3,1))
symm1 = [H; symm zeros(size(symm)); zeros(size(symm)) symm];
Z = null(full(symm1));
Temp = [reshape(TOR(:,:,:,1),N^3,1); reshape(POL(:,:,:,1),N^3,1)];
Temp = Z*(Z\Temp);
TOR(:,:,:,1) = reshape(Temp(1:N^3),N,N,N);
POL(:,:,:,1) = reshape(Temp(N^3+1:2*N^3),N,N,N);
% Kronecker up matrices and resize to N^2-2N+3
LHSkron1 = kron(C02(1:N-2,:),A(1:N-2,:))+kron(D2(1:N-2,:),coF(1:N-2,:));
C02big = kron(C02(1:N-2,:),C02(1:N-2,:));
Rbig = kron(C02(1:N-2,:),coF(1:N-2,:));

LHSkron1 = kron(speye(N),Rbig-alph*LHSkron1)-alph*kron(Dt^2,C02big);
LHSkron = kron(speye(2),LHSkron1);
%   4*N^2 - 6*N + 4 =  N*(4*N-6)+4 boundary rows.
% We need 2*(4N^2 - 4N) boundary rows in total.
symm = .5*kron(spdiags((1-(-1).^(1:N)).',0,N,N),kron(speye(N),(speye(N)-flip(speye(N)))*iDCT));
symm = symm + .5*kron(spdiags((1+(-1).^(1:N)).',0,N,N),kron(speye(N),(speye(N)+flip(speye(N)))*iDCT));
symm = [symm; kron([ones(N-1,1) -speye(N-1)]*iDFT,kron(speye(N),[zeros(1,(N-1)/2) 1 zeros(1,(N-1)/2)]*iDCT))];
symm1 = [symm zeros(size(symm)); zeros(size(symm)) symm];

[LHSkron, Indices] = RemoveDependencies([LHSkron; H; symm1]);
idx = find(Indices <= (N-2)^2*2*N,1,'last');
Indices = Indices(1:idx);
LHSkron = [LHSkron; zeros(2*N^3-size(LHSkron,1),2*N^3)];
LHSkron = LHSkron(1:2*N^3,:);
% column-pivoted QR factorization, 
%     LHSkron*E = Q*R 
[Q, R, E] = qr(full(LHSkron));
idx = find(abs(diag(R))>1e-14,1,'last');
for d = 1:T
    RHStor = zeros(N,N,N);
    RHSpol = zeros(N,N,N);
    for kk = 1:N
        RHStor(:,:,kk) = coF*TOR(:,:,kk,d)*C02.';
        RHSpol(:,:,kk) = coF*POL(:,:,kk,d)*C02.';
    end
    RHStor = RHStor(1:N-2,1:N-2,:);
    RHSpol = RHSpol(1:N-2,1:N-2,:);
    
    RHS = [reshape(RHStor,N*(N-2)^2,1); reshape(RHSpol,N*(N-2)^2,1)];
    RHS = RHS(Indices);
    RHS = [RHS; zeros(2*N^3-size(RHS,1),1)];%zeros(4*N^2 - 6*N + 4,1); zeros(extra,1)];
    
    % LHSkron*X = RHS
    % (LHSkron*E)*E'*X = RHS
    % (Q*R)*E'*X = RHS,  Q'*Q = I 
    % R*Y = Q'*RHS,  Y = E'*X  =>  X = E*Y. 
    
    Y = zeros(size(Q,1),1);
    Y(1:idx) = R(1:idx,1:idx) \ (Q(:,1:idx)'*RHS);
    X = E*Y;
    KapproxTOR = X(1:N^3);  % [LHSkron1; boundtest2]\[reshape(RHStor,N*(N-2)^2,1); zeros(N^3-N*(N-2)^2,1)];
    KapproxPOL = X(N^3+1:2*N^3);
    TOR(:,:,:,d+1) = reshape(KapproxTOR,N,N,N);
    POL(:,:,:,d+1) = reshape(KapproxPOL,N,N,N);
    %%%%%%%%%%%%%%%%%%%%%%%% 
end

%%
TORval=zeros(N,N,N,d);
POLval=zeros(N,N,N,d);
for d = 1:T
    TORval(:,:,:,d) = C2V_cyl(TOR(:,:,:,d),'rzt');
    POLval(:,:,:,d) = C2V_cyl(POL(:,:,:,d),'rzt');
end


%% Slice plot of the results

TimePlot_cyl(POLval,.1,alph,.0001);



%%
function [C1, E] = zeroDOF(C1, C2, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    for kk = 1:size(C1, 1)
        if ( abs(C1(kk,ii)) > 10*eps )
            c = C1(kk, ii); % Constant required to zero entry out.
            C1(kk,:) = C1(kk,:) - c*B(ii,:);
            E(kk,:) = E(kk,:) - c*G(ii,:)*C2.';
        end
    end
end

end
