
% A fast Poisson solver for the cylinder with zero homogeneous
% Dirichlet conditions.
%
% Take cylinder to be on (r,th,z) \in [0,1]x[-pi,pi]x[-1,1].
% We double up (x1) to get
%  (r,th,z) \in [-1,1]x[-pi,pi]x[-1,1].

% Discretization size:
n = 41;   % select odd

% Exact solution
k0 = besselroots(2,1);
[rr, tt, zz] = ndgrid( chebpts(n,[-1,1]), trigpts(n,[-pi,pi]), chebpts(n,[-1,1]) );
exact = @(r,th,z) besselj(2,k0*r).*sin(pi*z);
EXACT = exact(rr, tt, zz);
% Convert to Chebyshev--Fourier--Chebyshev coefficients:
for j = 1 : n
    EXACT(:,:,j) = chebtech2.vals2coeffs( EXACT(:,:,j) );
    EXACT(:,:,j) = trigtech.vals2coeffs( EXACT(:,:,j).' ).';
end
for j = 1 : n
    vj = reshape( EXACT( :, j, :), n, n );
    vj = chebtech2.vals2coeffs( vj.' ).';
    EXACT( :, j, : ) = reshape( vj, n, 1, n );
end

%  Right-hand side:
F = -EXACT;

% Discretize
% r^2*u_rr + r*u_r + u_{th,th} + r^2u_zz = r^2*f
Y = zeros(n, n, n);
X = zeros(n, n, n);

M2 = ultraS.multmat( n, [.5;0;-.5], 2 );
Mr = ultraS.multmat( n, [0;1], 1 );
D1 = ultraS.diffmat( n, 1);
D2 = ultraS.diffmat( n, 2);
S01 = ultraS.convertmat( n, 0, 1);
S11 = ultraS.convertmat( n, 1, 1);
B = S01 \ ( M2*D2 - 4*S11*Mr*D1 - 2*S01 );
% Construct useful matrices: 
j = (0:n-1)';
dsub = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j+2);
dsup = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j);
d = -dsub - dsup;
M = spdiags([dsub d dsup], [-2 0 2], n, n);   % Multiplication matrix for (1-x^2) in the C^(3/2) basis
D_inv = spdiags(-1./(j.*(j+3)+2), 0, n, n);     % D^{-1} undoes the scaling from the Laplacian identity
T = D_inv * M;                                 % Construct T = D^{-1} * M 
S = Chebyshev2ultra(n);
B = S \ T; 

M1 = ultraS.multmat( n, [.125;0;0;0;-.125], 2 );
M2 = ultraS.multmat( n, [.5;0;-.5], 0 );
M3 = ultraS.multmat( n, [0;.25;0;-.25], 1 );
Mr3 = ultraS.multmat( n, [0;.75;0;.25], 1 );
Mr1 = ultraS.multmat( n, [0;1], 1 );
Mr2 = ultraS.multmat( n, [.5;0;.5], 0 );
D1 = ultraS.diffmat( n, 1 );
D2 = ultraS.diffmat( n, 2 );
S01 = ultraS.convertmat( n, 0, 1 );
S11 = ultraS.convertmat( n, 1, 1 );
C = S01 \ (M1*D2 - 4*S11*Mr3*D1 - 2*S01*Mr2 + S11*M3*D1 - 2*S01*Mr2 ); 

S = leg2cheb( ultra1mx2Legendre( n ) );

eA = eig(full(M2\C));
a = max(eA);
b = min(eA);
eB = eig(full(B));
c = max(eB);
d = min(eB);
tic
% Fourier diff is diagonal, so 3D PDE separates into n decoupled 2D PDEs:
for k = -(n-1)/2 : (n-1)/2
    
    % Solve r^2*w_rr + r*w_r + r^2w_zz - k^2*w = r^2*f:
    
    % Build A:
    A = C - k^2 * M2;
    ak = a - k^2; 
    bk = b - k^2;
    if (k == 0 )
       keyboard 
    end
    
    Mr2 = ultraS.multmat( n, [.5;0;.5], 0 );
    Fk = reshape(F(:, k+(n+1)/2, :),n,n);
    Yk = lyap(M2\A, B.', -Fk );   % O(n^3)
    Xk = Sylvester_adi( M2\A, -B.', Fk, ak, bk, c, d); % This is O(n^2 log(1/eps))
    
    norm( Xk - Yk ) 
    Yk = S * Yk; 
    Xk = S * Xk; 
    Y(:, k+(n+1)/2, :) = reshape( Yk, n, 1, n );
    X(:, k+(n+1)/2, :) = reshape( Xk, n, 1, n );
end
toc

r = chebfun( 'x' );
rr = 1-r.^2;
M = ultraS.multmat( n, rr, 0 );
for k = 1:n
    err = M*EXACT(:,:,1) - Y(:,:,1);
    myerror(k) = norm(err(:));
end
norm(myerror)
