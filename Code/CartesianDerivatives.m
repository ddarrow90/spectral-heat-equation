N=16;
Ops = Operators(N,'D0 R0');
D = Ops.D0;
R = Ops.R0;

f = @(r,z,t) r.^2;
FF = map2grid_cyl(f,N);

%% D_x

GG = zeros(N,N,N+1);
for j = -N/2+1:N/2-1
    GG(:,:,j+N/2+1)=.5*D*FF(:,:,j+N/2)+.5*D*FF(:,:,j+N/2+2) ... % cos(th)D_r
        - R\(.5*(j-1)*FF(:,:,j+N/2)-.5*(j+1)*FF(:,:,j+N/2+2)); % - sin(th)D_th/R
end
GG(:,:,1) = .5*D*FF(:,:,2) + R\(.5*(-N/2+1)*FF(:,:,2));
GG(:,:,N+1) = .5*D*FF(:,:,N) - R\(.5*(N/2-1)*FF(:,:,N));

%% Horizontal Laplacian

for j = -N/2+1:N/2-1
    GG(:,:,j+N/2+1)=R\(D*R*D*FF(:,:,j+N/2+1))+(R^2)\(-j^2*FF(:,:,j+N/2+1));
end

%% Times sin(th)

N=16;
f = @(r,z,t) r.^2;
FF = func2grid(f,N,'rzt',1);
GG = zeros(N,N,N);
for j = -(N-1)/2+1:(N-1)/2-1
    GG(:,:,j+(N+1)/2)=-.5*1i*FF(:,:,j-1+(N+1)/2)+.5*1i*FF(:,:,j+1+(N+1)/2);
end

GG(:,:,1) = .5*1i*FF(:,:,2);
GG(:,:,N) = -.5*1i*FF(:,:,N-1);