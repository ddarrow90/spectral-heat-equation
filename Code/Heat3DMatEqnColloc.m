function CFS = Heat3DMatEqnColloc(f,N,T,alph)
% A fast heat equation solver for the cylinder.
%
% Take cylinder to be on (th,r,z,t) in [-pi,pi]x[0,1]x[-1,1]x[0,1].
% 
% We double up in the "r"-variable to get
% 
%  (th, r, z) in [-pi,pi]x[-1,1]x[-1,1]
%
% and use a Fourier-Chebyshev-Chebyshev expansion. 
% 
% Author: David Darrow, February 2017.


% Cylindrical harmonic
% f = @(t,r,z) feval(diskfun.harmonic(4,3),t,r,'polar').*sin(2*pi*z); %r.^3.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

% Construct grid (t,r,z) = Fourier-Chebyshev-Chebyshev:
t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );
[rr, zz, tt] = ndgrid(r, z, t);


Ops = Operators(N,'C02, D2, A, coF, DCT, iDCT');
DCT = Ops.DCT;
iDCT = Ops.iDCT;
C02 = full(iDCT*Ops.C02*DCT);
D2 = full(iDCT*Ops.D2*DCT);
A = full(iDCT*Ops.A*DCT);
coF = full(iDCT*Ops.coF*DCT);

% Store values of initial function in (p, r, z) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
CFS = zeros(size(zz));
for j=1:N
    CFS(:,j,:) = f(rr(:,j,:),zz(:,j,:), tt(:,j,:));
end

% Don't actually convert f from value space to coefficient space in
% Chebyshev. Only in Fourier.
CFS(:,:,:,1) = V2C_cyl(CFS,'00t');
CFS(:,:,:,2:T)=zeros;

Hist = zeros(N,N,N,30);

%1/10 time step; BDF2 30 times!


RHS = CFS(:,:,:,1);
q = alph/10;
Hist(:,:,:,1) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);%Heat3DSolver(cfs,Q,N,C02,A,D2,coF)


for p = 1:29
        RHS = Hist(:,:,:,p);
        Hist(:,:,:,p+1) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

CFS(:,:,:,[2 3 4])=Hist(:,:,:,[10 20 30]);

%t=5:T; full time step
q = (12/25)*alph;
for d = 5:T
        RHS = (48/25)*CFS(:,:,:,d-1)-(36/25)*CFS(:,:,:,d-2)+(16/25)*CFS(:,:,:,d-3)-(3/25)*CFS(:,:,:,d-4);
        CFS(:,:,:,d) = Heat3DSolver(RHS,q,N,C02,A,D2,coF);
end

% Convert back to value space
CFS = C2V_cyl(CFS,'00t0');


end

function Fapprox = Heat3DSolver(cfs,Q,N,C02,A,D2,coF)

% Solves the screened Poisson equation, [I-Q*Del^2]u=CFS in coeff. space
% Homogeneous Dirichlet boundary conditions.

% Discretization size: N

Fapprox = zeros(N,N,N);
Ops = Operators(N,'DCT');
DCT = Ops.DCT;

bc = .05*[ (-1).^(0:N-1) > 0 ; 
       (-1).^(0:N-1) < 0 ]*DCT; 

for k = -(N-1)/2 : (N-1)/2
    
    AA = coF - Q*(A - k^2*C02);
    BB = C02; 
    CC = -Q*coF; 
    DD = D2; 
    RHS = coF*cfs(:,:,k+(N+1)/2)*C02.';
    %RHS = RHS(1:N-2, 1:N-2);
    
    [AA, RHS] = zeroDOF(AA, BB, RHS, bc, zeros(2,N));
    [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, zeros(2,N));
    RHS = RHS.';
    [CC, RHS] = zeroDOF(CC, DD, RHS, bc, zeros(2,N));
    [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, zeros(2,N));
    RHS = RHS.';
    
    AA = AA(2:end-1,2:end-1);
    BB = BB(2:end-1,2:end-1);
    CC = CC(2:end-1,2:end-1);
    DD = DD(2:end-1,2:end-1);
    RHS = RHS(2:end-1,2:end-1);
    
    X22 = bartelsStewart(AA, BB, CC, DD, RHS );
    X = [ zeros(1,N) ; zeros(N-2,1) X22 zeros(N-2,1); zeros(1,N)]; 
    Fapprox(:, :, k+(N+1)/2) = X;

end

end

function [C1, E] = zeroDOF(C1, C2, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.
sampleends = [1 size(C1, 2)];
for ii = 1:2 % For each boundary condition, zero a column.
    for kk = 1:size(C1, 1)
        if ( abs(C1(kk,ii)) > 10*eps )
            c = C1(kk, sampleends(ii)); % Constant required to zero entry out.
            C1(kk,:) = C1(kk,:) - c*B(ii,:);
            E(kk,:) = E(kk,:) - c*G(ii,:)*C2.';
        end
    end
end

end
