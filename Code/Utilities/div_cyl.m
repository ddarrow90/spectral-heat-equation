function X = div_cyl(Fr,Ft,Fz)
% Given R, Theta, and Z components of a cylindrical vector field
N = size(Fr,1);
Ops = Operators(size(Fr,1),'D0 R0');
D0 = Ops.D0;
R0 = Ops.R0;

X = zeros(N,N,N+1);
for k = -N/2 : N/2
    X(:,:,k+N/2+1) = D0*Fr(:,:,k+N/2+1) + R0\(Fr(:,:,k+N/2+1));
    X(:,:,k+N/2+1) = X(:,:,k+N/2+1)+1i*k*(R0\(Ft(:,:,k+N/2+1)));
    X(:,:,k+N/2+1) = X(:,:,k+N/2+1)+Fz(:,:,k+N/2+1)*D0.';
end
end