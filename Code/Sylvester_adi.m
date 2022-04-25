function X = Sylvester_adi( A, B, F, a, b, c, d)
% Solve A*X-X*B = F using ADI. 
%
% eig(A) \in [a,b] and eig(B) \in [c,d]. 
% 
% Example: 
%
% n = 40;
% U = qr(rand(n));
% A = U * diag(2*rand(n,1)-1) / U;
% U = qr(rand(n));
% B = U * diag(rand(n,1)+10) / U;
% F = rand(n);
% 
% eA = eig( A );
% a = min( eA );
% b = max( eA );
% 
% eB = eig( B );
% c = min( eB );
% d = max( eB );
% 
% Yk = lyap(A, -B, -F );
% X = Sylvester_adi( A, B, F, a, b, c, d);
% norm( X - Yk )
% 

m = size(A,1);
n = size(B,1);

bb = abs(c-a).*abs(d-b)./abs(c-b)./abs(d-a);
gam = -1 + 2*bb + 2*sqrt(bb^2-bb);
z1 = -gam; z2 = -1; z3 = 1;   % [-gam, -1] \cup [1, gam] 
w1 = b; w2 = a; w3 = c;       % [a, b] \cup [c, d]  
AA = det( [ z1*w1 w1 1 ; z2*w2 w2 1 ; z3*w3 w3 1] );
BB = det( [ z1*w1 z1 w1 ; z2*w2 z2 w2 ; z3*w3 z3 w3] );
CC = det( [z1 w1 1 ; z2 w2 1 ; z3 w3 1] );
DD = det( [ z1*w1 z1 1 ; z2*w2 z2 1 ; z3*w3 z3 1] );
T = @(z) (AA*z+BB)./(CC*z+DD);

tol = 1e-13;
cos2beta = 2/(1+.5*(1/gam+gam/1));
M = 2/cos2beta - 1;
kprime = 1/(M+sqrt(M^2-1));
[Kprime, K] = ellipk(kprime);
J = ceil(K/2/Kprime/pi*log(4/tol) );
[~, ~, dn] = ellipj((1:2:(2*J-1))*K/2/J, sqrt(1-kprime.^2));
p = sqrt(1*gam/kprime)*dn;                  % OPTIMAL SHIFT PARAMETERS!
q = T(p);
p = T(-p);

X = zeros(m, n); In = speye(n); Im = speye(m);
for j = 1:numel( p )
    X = (A-q(j)*In) \ ( F + X*(B-q(j)*Im) );
    X = ( (A-p(j)*In)*X - F ) / (B-p(j)*Im);
end

end