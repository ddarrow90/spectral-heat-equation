function [Fr,Ft,Fz] = reversePT2(T,P)
N = size(P,1);

Ops = Operators(N,'D0,R0');
Fr = zeros(N,N,N);
Ft = zeros(N,N,N);
Fz = zeros(N,N,N);

for j = -(N-1)/2:(N-1)/2
    R_temp = 1i*j*P(:,:,j+(N+1)/2); % Multiplied by R here.
    Th_temp = -Ops.D0*P(:,:,j+(N+1)/2);
    
    Fr(:,:,j+(N+1)/2) = 1i*j*Ops.R0*T(:,:,j+(N+1)/2) - Ops.R0^2*Th_temp*Ops.D0.'; % And here.
    Ft(:,:,j+(N+1)/2) = - Ops.R0^2*Ops.D0*T(:,:,j+(N+1)/2) + Ops.R0*R_temp*Ops.D0.'; % Already half done here.
    Fz(:,:,j+(N+1)/2) = Ops.R0^2*Ops.D0*Th_temp + Ops.R0*Th_temp - 1i*j*R_temp; % Once more here.
end
end