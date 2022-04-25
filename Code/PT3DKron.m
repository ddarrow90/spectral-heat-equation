function [Tor,Pol,Sol] = PT3DKron(CFSr,CFSt,CFSz,N,checkparam)

% Creates a PT decomposition for divergence free fields, such that 
% F = curl(Tor*e_3) + curl curl(Pol*e_3)

% Toroidal Scalar: -hLap(Tor)=e_3 * curl(F)
% Poloidal Scalar: -hLap(Tor)=e_3 * F
% Where e_3 is the vertical unit vector and hLap is the Laplacian on the
% disk.

% Checkparam indicates whether or not to check for solenoidal-ness. "Sol"
% indicates this, being 1 for incompressible fields and 0 otherwise.

%{
Fr = @(t,r,z) 4.*r.*z+0.*t;
Ft = @(t,r,z) 2.*r.*(1-z.^2)+0.*t;
Fz = @(t,r,z) 4.*(1-z.^2)+0.*t;
N = 15;
%}
if ( nargin < 5 )
    checkparam = 0;
end

% Make sure the parity of N is even:
if ( mod(N, 2) )  
    N = N + 1;
end

Ops = Operators(N,'C01 C12 C02 D1 D0 R0 R2 R1 A coF');
C01 = Ops.C01;
C12 = Ops.C12;
C02 = Ops.C02;
D1 = Ops.D1;
R2 = Ops.R2;
R1 = Ops.R1;
R0 = Ops.R0;
D0 = Ops.D0;
A = Ops.A;
coF = Ops.coF;

if checkparam
    if (infnorm_cyl(div_cyl(Fr,Ft,Fz)) > N*eps)
        Sol = 0;
        max(real(check))
    end
end
    
Tor = zeros(N,N,N+1);
Pol = zeros(N,N,N+1);
bc = [ (-1).^(0:N-1) > 0 ; 
       (-1).^(0:N-1) < 0 ]; 
%bc1 = bc/C01;

for k = -N/2 : N/2
    RHS = sparse(N^2,N^2);
    AAA = sparse(N^2,N^2);
    
    AA = A - k^2*C02; % Good.
    RHSt = -(R2*C02+R2^2*C12*D1)*CFSt(:,:,k+N/2+1)+k*1i*R2*C02*CFSr(:,:,k+N/2+1); % Good.
    RHSp = -coF*CFSz(:,:,k+N/2+1); % Good.
    AA = AA(1:end-2,:);
    AA = kron(speye(N),AA);
    RHSt = RHSt(1:N-2,:);
    RHSt = reshape(RHSt,N^2-2*N,1);
    RHSp = RHSp(1:N-2,:);
    RHSp = reshape(RHSp,N^2-2*N,1); % Good up to here.
    %[AA, RHS] = zeroDOF(AA, RHS, bc, zeros(2,N));
    
    bAAA = [1i*k*kron(speye(N),bc) kron(D0,bc*R0*D0)];%kron(speye(N),bc1)*
    bRHS = bc*R0*CFSr(:,:,k+N/2+1);
    bRHS = reshape(bRHS,2*N,1); % Good up to here.
    %bRHS = zeros(2*N,1);
    
    bAAA1 = [-kron(speye(N),bc*R0*D0) 1i*k*kron(D0,bc)]; % Boundary for each thing?
    bRHS1 = bc*R0*CFSt(:,:,k+N/2+1);
    
    bRHS1 = reshape(bRHS1,2*N,1); % Good up to here.
    
    RHS = [RHSt; RHSp; bRHS; bRHS1];
    AAA = [kron(speye(2),AA); bAAA; bAAA1];
    %AAA = AAA(:,1:end-2*N-1);
    [Q, R, E] = qr(full(AAA));
    idx = find(abs(diag(R))>1e-14,1,'last');
    
    % AAA*E*E'*X = RHS
    % R*E'*X = Q'*RHS
    Y = zeros(size(Q,1),1);
    Y(1:idx) = R(1:idx,1:idx) \ (Q(:,1:idx)'*RHS);
    X = E*Y; 
    
%     X = AAA\RHS;
    XTor = X(1:N^2);
    XPol = X(N^2+1:2*N^2);
    %spy(X)
    Tor(:, :, k+N/2+1) = reshape(XTor,N,N); %FapproxTor
    Pol(:, :, k+N/2+1) = reshape(XPol,N,N); %FapproxPol
end
end

% for k = -(N-1)/2 : (N-1)/2
% infnorm_cyl(i*k*bc*T(:,:,k+N/2+1) + bc*R0*D0*P(:,:,k+N/2+1)*D0.'-bc*R0*CFSr(:,:,k+N/2+1))
% end
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