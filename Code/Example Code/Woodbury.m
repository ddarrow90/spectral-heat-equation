% Test problem: 
n = 10; k = 2; 
A = rand(n,n);
b = rand(n,1); 

% Low rank update of A: 
U = rand(n,k); 
V = rand(k,n); 
I = eye(k,k); 
B = A + U*V; 

% We want to solve, Bx = b by Woodbury: 
y = A \ b; 
x = y - A\((U*((I + V*(A\U))\V))*y);

% Check: 
exact = B \ b; 
norm( exact - x )