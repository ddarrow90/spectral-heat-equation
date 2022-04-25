function FFapprox = ScrPoisson3DMatEqn(f,N,alph)

% Solves the screened Poisson equation, [Del^2-1/alph]u=-f/alph
% Or, identically, [I-alph*Del^2]u=f

% Discretization size: N

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
FF = f(tt,rr,zz);

% Convert f from value space to coefficient space
CFS = reshape(FF,N,N^2);
CFS = chebtech2.vals2coeffs(CFS);
CFS = reshape(CFS,N,N,N);
CFS = permute(CFS,[3 2 1]);
CFS = reshape(CFS,N,N^2);
CFS = chebtech2.vals2coeffs(CFS);
CFS = reshape(CFS,N,N,N);
CFS = permute(CFS,[2 3 1]);
CFS = reshape(CFS,N,N^2);
CFS = trigtech.vals2coeffs(CFS);
CFS = reshape(CFS,N,N,N);
CFS = permute(CFS,[2 3 1]);

Fapprox = zeros(N,N,N);
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
    
    AA = coF - alph*(A - k^2*U02);
    BB = U02; 
    CC = -alph*coF; 
    DD = D2; 
    RHS = coF*CFS(:,:,k+(N+1)/2)*U02.';
    %RHS = RHS(1:N-2, 1:N-2);
    
    [AA, RHS] = zeroDOF(AA, BB, RHS, bc, zeros(2,N));
    [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, zeros(2,N));
    RHS = RHS.';
    [CC, RHS] = zeroDOF(CC, DD, RHS, bc, zeros(2,N));
    [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, zeros(2,N));
    RHS = RHS.';
    
    AA = AA(1:end-2,3:end);
    BB = BB(1:end-2,3:end);
    CC = CC(1:end-2,3:end);
    DD = DD(1:end-2,3:end);
    RHS = RHS(1:end-2,1:end-2);
    
    X22 = bartelsStewart(AA, BB, CC, DD, RHS );
    X12 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X22);
    X21 = ( bc(1:2,1:2).' \ (-bc(1:2,3:end)*X22.') ).';
    X11 = bc(1:2,1:2) \ (-bc(1:2,3:end)*X21);
    X = [ X11 X12 ; X21 X22 ]; 
    Fapprox(:, :, k+(N+1)/2) = X;

end

FFapprox = reshape(Fapprox,N,N^2);
FFapprox = chebtech2.coeffs2vals(FFapprox);
FFapprox = reshape(FFapprox,N,N,N);
FFapprox = permute(FFapprox,[2 1 3]);
FFapprox = reshape(FFapprox,N,N^2);
FFapprox = chebtech2.coeffs2vals(FFapprox);
FFapprox = reshape(FFapprox,N,N,N);
FFapprox = permute(FFapprox,[3 2 1]);
FFapprox = reshape(FFapprox,N,N^2);
FFapprox = real(trigtech.coeffs2vals(FFapprox));
FFapprox = reshape(FFapprox,N,N,N);
FFapprox = permute(FFapprox,[2 1 3]);

end

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