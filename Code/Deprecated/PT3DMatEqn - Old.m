function [Tor,Pol,Sol] = PT3DMatEqn(Fr,Ft,Fz,N,checkparam)

% Creates a PT decomposition for divergence free fields, such that 
% F = curl(Tor*e_3) + curl curl(Pol*e_3)

% Toroidal Scalar: -hLap(Tor)=e_3 * curl(F)
% Poloidal Scalar: -hLap(Tor)=e_3 * F
% Where e_3 is the vertical unit vector and hLap is the Laplacian on the
% disk.

% Checkparam indicates whether or not to check for solenoidal-ness. "Sol"
% indicates this, being 1 for incompressible fields and 0 otherwise.

% Test: F = 4rzR + 2r(1-z^2)T+4(1-z^2)Z
% Answer: Tor = Pol = (1-r^2)(1-z^2)
% Discretization size: N

%Fr = @(t,r,z) 4.*r.*z+0.*t;
%Ft = @(t,r,z) 2.*r.*(1-z.^2)+0.*t;
%Fz = @(t,r,z) 4.*(1-z.^2)+0.*t;
%N = 50;
if ( nargin < 5 )
    checkparam = 0;
end

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

% Construct grid (t,r,z) = Fourier-Chebyshev-Chebyshev:
t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );
[tt, rr, zz] = meshgrid(t, r, z);

% Construct useful ultraspherical operators: 
% Conversion operator: (ChebT -> ChebU)
U01 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); 
U01(1, 1) = 1; 

% Conversion operator: (ChebU -> C^2)
K = 1./(1:N)';
U12 = spdiags([K -K], [0 2], N, N);

% Conversion operator: (ChebT -> C^2)
U02 = U12*U01;

% First-order diff: (ChebT -> ChebU)
D1 = spdiags((0:N)', 1, N, N); 

% Second-order diff: (ChebT -> C^2)
D2 = spdiags(2*(0:N)', 2, N, N); 

% Multiplication by "r": (in C^2)
K = (1:N)'./(4:2:2*N+2)';
K1 = (3:N+2)'./(4:2:2*N+2)';
R = spdiags([K K1], [-1 1], N, N); 

% Construction of "r"-part of Laplacian on cylinder:
A = R^2*D2 + R*U12*D1;   % A[u] = r^2*u_rr + r*u_r

% Construction of "f"-part of Laplacian on cylinder
coF=R^2*U02;

% Store values of test functions in (t, r, z) grid
FFr = Fr(tt,rr,zz);
FFt = Ft(tt,rr,zz);
FFz = Fz(tt,rr,zz);

% Convert f from value space to coefficient space
CFSz = reshape(FFz,N,N^2);
CFSz = chebtech2.vals2coeffs(CFSz);
CFSz = reshape(CFSz,N,N,N);
CFSz = permute(CFSz,[3 2 1]);
CFSz = reshape(CFSz,N,N^2);
CFSz = chebtech2.vals2coeffs(CFSz);
CFSz = reshape(CFSz,N,N,N);
CFSz = permute(CFSz,[2 3 1]);
CFSz = reshape(CFSz,N,N^2);
CFSz = trigtech.vals2coeffs(CFSz);
CFSz = reshape(CFSz,N,N,N);
CFSz = permute(CFSz,[2 3 1]);

CFSr = reshape(FFr,N,N^2);
CFSr = chebtech2.vals2coeffs(CFSr);
CFSr = reshape(CFSr,N,N,N);
CFSr = permute(CFSr,[3 2 1]);
CFSr = reshape(CFSr,N,N^2);
CFSr = chebtech2.vals2coeffs(CFSr);
CFSr = reshape(CFSr,N,N,N);
CFSr = permute(CFSr,[2 3 1]);
CFSr = reshape(CFSr,N,N^2);
CFSr = trigtech.vals2coeffs(CFSr);
CFSr = reshape(CFSr,N,N,N);
CFSr = permute(CFSr,[2 3 1]);

CFSt = reshape(FFt,N,N^2);
CFSt = chebtech2.vals2coeffs(CFSt);
CFSt = reshape(CFSt,N,N,N);
CFSt = permute(CFSt,[3 2 1]);
CFSt = reshape(CFSt,N,N^2);
CFSt = chebtech2.vals2coeffs(CFSt);
CFSt = reshape(CFSt,N,N,N);
CFSt = permute(CFSt,[2 3 1]);
CFSt = reshape(CFSt,N,N^2);
CFSt = trigtech.vals2coeffs(CFSt);
CFSt = reshape(CFSt,N,N,N);
CFSt = permute(CFSt,[2 3 1]);

if checkparam
    r2Div = zeros(N,N,N);
    Sol = 1;
    for k = -(N-1)/2 : (N-1)/2
        r2Div(:,:,k+(N+1)/2) = (R*U02+R^2*U12*D1)*CFSr(:,:,k+(N+1)/2)*U02.';
        r2Div(:,:,k+(N+1)/2) = r2Div(:,:,k+(N+1)/2)-k^2*R^2*U02*CFSt(:,:,k+(N+1)/2)*U02.';
        r2Div(:,:,k+(N+1)/2) = r2Div(:,:,k+(N+1)/2)+R^2*U02*CFSz(:,:,k+(N+1)/2)*D1.'*U12.';
        check = reshape(r2Div(:,:,k+(N+1)/2),N^2,1);
        if (max(check) > 10*eps)
            Sol = 0;
        end
    end
end
    
FapproxTor = zeros(N,N,N);
FapproxPol = zeros(N,N,N);
% Kapprox = zeros(N^2,1);
%U02 = U02(1:N-2,:);
%A = A(1:N-2,:);
%D2 = D2(1:N-2,:);
%coF = coF(1:N-2,:);
% LHSkron1 = kron(U02, A) + kron(D2, coF);
% U02big = kron(U02, U02);
bc = [ (-1).^(0:N-1) > 0 ; 
       (-1).^(0:N-1) < 0 ]; 
% BCs = [kron(ones(1,N), eye(N)) ;
%      kron((-1).^(1:N), eye(N)) ; 
%      kron( eye(N-2,N), ones(1,N) ) ;
%      kron( eye(N-2,N), (-1).^(1:N) ) ];
for k = -(N-1)/2 : (N-1)/2
    
    AA = A - k^2*U02;
    RHS = -(R^2*U12*D1+R*U02)*CFSt(:,:,k+(N+1)/2)+k*1i*R*U02*CFSr(:,:,k+(N+1)/2);
    %RHS = RHS(1:N-2, 1:N-2);
    
    [AA, RHS] = zeroDOF(AA, RHS, bc, zeros(2,N));
    
    AA = AA(1:end-2,3:end);
    RHS = RHS(1:end-2,1:end);
    
    X22 = AA\RHS;
    X12 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X22);
    X = [ X12 ; X22 ]; 
    FapproxTor(:, :, k+(N+1)/2) = X;
end

for k = -(N-1)/2 : (N-1)/2
    
    AA = A - k^2*U02;
    RHS = -coF*CFSz(:,:,k+(N+1)/2);
    %RHS = RHS(1:N-2, 1:N-2);
    
    [AA, RHS] = zeroDOF(AA, RHS, bc, zeros(2,N));
    %{
    AA = AA(1:end-2,3:end);
    RHS = RHS(1:end-2,1:end-2);
    
    X22 = AA\RHS;
    X12 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X22);
    X21 = ( bc(1:2,1:2).' \ (-bc(1:2,3:end)*X22.') ).';
    X11 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X21);
    X = [ X11 X12 ; X21 X22 ]; 
    %}
    AA = AA(1:end-2,3:end);
    RHS = RHS(1:end-2,1:end);
    
    X22 = AA\RHS;
    X12 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X22);
    X = [ X12 ; X22 ]; 
    FapproxPol(:, :, k+(N+1)/2) = X;

end

Pol = reshape(FapproxPol,N,N^2);
Pol = chebtech2.coeffs2vals(Pol);
Pol = reshape(Pol,N,N,N);
Pol = permute(Pol,[2 1 3]);
Pol = reshape(Pol,N,N^2);
Pol = chebtech2.coeffs2vals(Pol);
Pol = reshape(Pol,N,N,N);
Pol = permute(Pol,[3 2 1]);
Pol = reshape(Pol,N,N^2);
Pol = real(trigtech.coeffs2vals(Pol));
Pol = reshape(Pol,N,N,N);
Pol = permute(Pol,[2 1 3]);

Tor = reshape(FapproxTor,N,N^2);
Tor = chebtech2.coeffs2vals(Tor);
Tor = reshape(Tor,N,N,N);
Tor = permute(Tor,[2 1 3]);
Tor = reshape(Tor,N,N^2);
Tor = chebtech2.coeffs2vals(Tor);
Tor = reshape(Tor,N,N,N);
Tor = permute(Tor,[3 2 1]);
Tor = reshape(Tor,N,N^2);
Tor = real(trigtech.coeffs2vals(Tor));
Tor = reshape(Tor,N,N,N);
Tor = permute(Tor,[2 1 3]);
end

function [C1, E] = zeroDOF(C1, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    for kk = 1:size(C1, 1)
        if ( abs(C1(kk,ii)) > 10*eps )
            c = C1(kk, ii); % Constant required to zero entry out.
            C1(kk,:) = C1(kk,:) - c*B(ii,:);
            E(kk,:) = E(kk,:) - c*G(ii,:);
        end
    end
end

end