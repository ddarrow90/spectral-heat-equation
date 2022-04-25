function [Fr,Ft,Fz] = reversePT(T,P)
N = size(P,1);

Ops = Operators(N,'D0,R0');
R0 = Ops.R0;
D0 = Ops.D0;
Fr = zeros(N,N,N+1);
Ft = zeros(N,N,N+1);
Fz = zeros(N,N,N+1);

for j = -N/2:N/2
    R_temp = R0\(1i*j*P(:,:,j+N/2+1));
    Th_temp = -D0*P(:,:,j+N/2+1);
    
    Fr(:,:,j+N/2+1) = R0\(1i*j*T(:,:,j+N/2+1)) - Th_temp*D0.';
    Ft(:,:,j+N/2+1) = - D0*T(:,:,j+N/2+1) + R_temp*D0.';
    Fz(:,:,j+N/2+1) = D0*Th_temp + R0\(Th_temp - 1i*j*R_temp);
end

end