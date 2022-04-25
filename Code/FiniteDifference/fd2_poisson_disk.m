%%
% Size of the grid; n should be even so that r=0 is not a grid point.
n = 512; n = n+mod(n,2); m = n;
r = linspace(-1,1,n+2)'; r = r(2:n+1); hr = 2/(n+1);
th = pi*trigpts(m);
[tt,rr] = meshgrid(th,r);

% Second order accurate differentiation matrices for d/dr and d^2/dr^2
D = (1/2/hr)*spdiags(ones(n,1)*[-1 1],[-1 1],n,n);
D2 = (1/hr/hr)*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
% Multiplication matrices in value space by R and R^2
R = diag(r); R2 = R.^2;
I = speye(n);

% Test problem: exact solution to recover (needs zero boundary conditions)
u_exact = diskfun(@(th,r) (1-r.^2).*exp(r.^2.*cos(th-0.1).*sin(th+0.3)),'polar');
% laplacian(u) = f
f = laplacian(u_exact);
% Sample f at the interior grid points
f = f(tt,rr,'polar');         

% Fourier transform in theta to decouple theta and r variables.
F = trigtech.vals2coeffs(f.').';
U = 0*F;
% Fourier modes to solve for, handle the zero mode separately
modes = -(m/2):(m/2-1); modes(m/2+1) = [];
for k = modes
    L = R2*D2 + R*D - k^2*I;  % Differentiation matrix for mode k.
    U(:,k + m/2+1) = L\(R2*F(:,k + m/2+1));
end
% Solve for the zero mode
U(:,m/2+1) = (R*D2 + D)\(R*F(:,m/2+1));
% Inverse Fourier transform in theta to get the solution.
u = real(trigtech.coeffs2vals(U.').');
% Exact solution on the grid
u_exact_grid = u_exact(tt,rr,'polar');

% Remove the solution corresponding to r < 0.
u = u(n/2+1:end,:); u_exact_grid = u_exact_grid(n/2+1:end,:);
tt = tt(n/2+1:end,:); rr = rr(n/2+1:end,:);
xx = rr.*cos(tt); yy = rr.*sin(tt);

% Fill in the theta = pi values so plotting on the disk looks correct:
u = [u u(:,1)]; u_exact_grid = [u_exact_grid u_exact_grid(:,1)];
xx = [xx xx(:,1)]; yy = [yy yy(:,1)];

% Compute and plot the error
err = u - u_exact_grid;
err_inf = norm(err,inf)/norm(u_exact_grid,inf)
surf(xx,yy,err), shading interp




