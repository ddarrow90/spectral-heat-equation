

%%
n = 200; 
tol = 1e-15;
m = n;
n = 2*n+1;

%get coeffs and set up RHS 
f = diskfun(@(x,y) cos(x.*y));

% Construct operators
D1 = ultraS.diffmat( n, 1 );              % 1st order ultraS diffmat
D2 = ultraS.diffmat( n, 2 );              % 2nd order ultraS diffmat
Mr = ultraS.multmat(n,[0;1],1);           % multiplication of r in ChebU 
Mr2 = ultraS.multmat(n,[.5;0;.5],2);      % multiplication of r^2 in ultra2
Mr2c = ultraS.multmat(n,[.5;0;.5],0);     % multiplication of r^2 in Cheb
S1 = ultraS.convertmat( n, 0, 1 );        % convert chebT coeffs -> ultra2
S12 = ultraS.convertmat( n, 1, 1);        % convert chebU coeffs -> ultra2

%left multiply operator
A = Mr2*(D2) + S12*Mr*D1;

%right multiply operator: B is D2 for trig
B = -1*((1:m/2).^2).';
B = [B(end:-1:1); 0; B(1:end-1)] ; 
B = spdiags(B, 0, m, m); 

% [u, ZZ, DD, YY] = poisson_disk_FIADI_final(f, n);
u = diskfun.poisson( f, @(th) 0*th, n, m ); 
EXACT = coeffs2( u, m, n ); 


F = S1*Mr2c*coeffs2( f, m, n); 
A = A - [zeros(n-2,n) ; A(end-1:end,:) ] + [zeros(n-2,n) ; ones(1,n) ; (-1).^(0:n-1) ];
% A(end,end) = A(end,end) - B(end,end); 
% A(end-1,end) = A(end-1,end-1) - B(end-1,end-1); 
F(end-1:end,:) = 0;

X = chebop2.bartelsStewart( A, eye(m,m), eye(n,n), B, F, 0, 0);
v = diskfun.coeffs2diskfun( X ); 


norm( X - EXACT)
 spy(abs(X-EXACT)>1e-13)
 
%% Disk Poisson solver, with partial regularity: 
n = 100; 
% N = chebop(@(r,u) r.^2.*diff((1-r.^2).*u, 2) + r.*diff((1-r.^2).*u));
N = chebop(@(r,u) diff((1-r.^2).*r.^2.*u, 2) + r.*(1-r.^2).*diff(u) + 2*(1-2*r.^2).*u);

A = zeros(n,n); 
M = ultraS.multmat( n, [.5;0;-.5], 0); 
S = ultraS.convertmat( n, 0, 0 ); 
for k = 1:n 
    v = chebpoly( k-1, [-1,1], 1 );
    A(:,k) = chebcoeffs( N(v), n, 'kind', 1);
    k
end
C = M \ A;

% check
[V, D] = eig( C ); 
cond( V )
plot(eig(C)+eps*1i, '.')

%% Run Poisson solver using ADI (right now just Bartel-Stewart): 
D = diag( (1i*(-n/2:n/2-1)).^2 ); 
f = diskfun(@(x,y) cos(x.*y));
F = coeffs2( f, n, n); 

u = diskfun.poisson( f, @(th) 0*th, n, n ); 
EXACT = coeffs2( u, n, n ); 

X = chebop2.bartelsStewart( A, eye(n,n), M, D, F, 0, 0);
Mr = ultraS.multmat( n, [.5;0;.5], 0);
X = Mr*M*X;
X(:,floor(n/2)+1) = EXACT(:,floor(n/2)+1);
v = diskfun.coeffs2diskfun( X ); 
norm( EXACT - X )

%% Poisson Cylinder, with partial regularity: 

n = 100; om = 3; 
N = chebop(@(r, u) r.^2.*(1-r.^2).*diff(u,2) - 5*r.^3.*diff(u) + r.^2.*om^2*u + 1*(-om^2*u + r.*diff(u))); 
N = chebop(@(r, u) r.^2.*(1-r.^2).*diff(u,2) + r.*(5-9*r.^2).*diff(u) + 4*u - 16.*r.^2.*u - om^2t(1-r.^2).*u );

A = zeros(n,n); 
M = ultraS.multmat( n, [.125;0;0;0;.125], 1 ); 
M = ultraS.multmat( n, [.5;0;-.5], 1);
S = ultraS.convertmat( n, 0, 0 ); 
for k = 1:n 
    v = chebpoly( k-1, [-1,1], 2 );
    A(:,k) = chebcoeffs( N(v), n, 'kind', 2);
    k
end
C = M \ A;
A = A(1:n,1:n); 
C = C(1:n,1:n); 

[V, D] = eig( C ); 
cond( V )
plot(eig(C)+eps*1i, '.')
% norm(C*C' - C'*C)

%%
spy(abs(A)>1e-3)

%%
n = 100;
N = chebop(@(x,u) diff((1-x.^2).*u,2));

A = zeros(n,n); 
for k = 1:n 
    v = chebpoly( k-1, [-1,1], 1 );
    A(:,k) = chebcoeffs( N(v), n, 'kind', 2);
    k
end

spy(abs(A)>1e-3)
[V, D] = eig( A ); 

% norm(C*C' - C'*C)


%%

g = trigtech( @(th) cos(8*pi*th) + sin(pi*th).^2); 
gc = (g.coeffs).';
gv = (g.values).'; 

n = length(g); 
X = zeros(n,n); 
X(1, floor(n/2)+1) = gc(1);

u = diskfun.coeffs2diskfun( X ); 

er = u(1,:) - chebfun(g,'periodic');
plot(er)


% ones(1,n)*C = g.values
% (-1).^(0:n-1)*C = fliplr(g.values)










