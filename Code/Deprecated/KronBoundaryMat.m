% Dave Darrow, 5/31/17.
%
% This code approximately multiplies two functions together, using a
% Chebyshev-Fourier-Chebyshev expansion and a Kronecker product.


f = @(t,r,z) r.^3.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
g = @(t,r,z) r.^3.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

% Discretization size, time discretization size: 
N = 20;

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

% Construct grid (t,r,z) = Fourier-Chebyshev-Chebyshev:
t = pi*trigpts( N );
r = chebpts( N );
z = chebpts( N );
[tt, rr, zz] = meshgrid(t, r, z);

%%
% Construct multiplication operator: 

% Store values of initial function in (p, r, z) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
FF = zeros(size(zz));
for j=1:N
    FF(:,:,j) = f(tt(:,:,j),rr(:,:,j),zz(:,:,j));
end

GG = zeros(size(zz));
for j=1:N
    GG(:,:,j) = f(tt(:,:,j),rr(:,:,j),zz(:,:,j));
end

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
FFc = reshape(CFS,N^2,N);
FFc = reshape(CFS,N^3,1);

CFS = reshape(GG,N,N^2);
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
GGc = reshape(CFS,N^2,N);
GGc = reshape(CFS,N^3,1);

DCT = chebtech2.coeffs2vals(eye(N));
DFT = trigtech.coeffs2vals(eye(N));


Mult = real(kron(DFT,kron(DCT,DCT))*GGc.*kron(DFT,kron(DCT,DCT))*FFc);

% Output1 = output in value space.
Output1 = reshape(Mult(:,:,:),N,N,N);
Output1 = permute(Output1,[3 1 2]);

%MultCoeff = output in coefficient space.
CFS = reshape(Output1,N,N^2);
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
MultCoeff = permute(CFS,[2 3 1]);
