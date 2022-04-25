% A fast heat equation solver for the cylinder.
%
% Take cylinder to be on (th,r,z,t) in [-pi,pi]x[0,1]x[-1,1]x[0,1].
% 
% We double up in the "r"-variable to get
% 
%  (th, r, z) in [-pi,pi]x[-1,1]x[-1,1]
%
% and use a Fourier-Chebyshev-Chebyshev expansion. 
% 
% Author: David Darrow, February 2017.

% Initial (time 0) Value: 
f = @(t,r,z) r.^3.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

% Cylindrical harmonic
%f = @(t,r,z) feval(diskfun.harmonic(4,3),t,r,'polar').*sin(2*pi*z);

% Discretization size, time discretization size: 
N = 9;
T = 15;
alph = .1;

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
% Construct useful ultraspherical operators: 
% Conversion operator: (ChebT -> ChebU)
C01 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); 
C01(1, 1) = 1; 

% Conversion operator: (ChebU -> C^2)
K = 1./(1:N)';
C12 = spdiags([K -K], [0 2], N, N);

% Conversion operator: (ChebT -> C^2)
C02 = C12*C01;

% First-order diff: (ChebT -> ChebU)
D1 = spdiags((0:N)', 1, N, N); 

% Second-order diff: (ChebT -> C^2)
D2 = spdiags(2*(0:N)', 2, N, N); 

% Multiplication by "r": (in C^2)
K = (1:N)'./(4:2:2*N+2)';
K1 = (3:N+2)'./(4:2:2*N+2)';
R = spdiags([K K1], [-1 1], N, N); 

% Construction of "r"-part of Laplacian on cylinder:
A = R^2*D2 + R*C12*D1;   % A[u] = r^2*u_rr + r*u_r

% Construction of "f"-part of Laplacian on cylinder
coF=R^2*C02;

% Store values of initial function in (p, r, z) grid
% Bug in diskfun.harmonic prevents this:
% FF = f(tt,rr,zz);
% Instead we have to loop over the zz levels
FF = zeros(size(zz));
for j=1:N
    FF(:,:,j) = f(tt(:,:,j),rr(:,:,j),zz(:,:,j));
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
CFS(:,:,:,2:T)=zeros;
Fapprox = zeros(N,N,N);
Kapprox = zeros(N^2,1);
C02 = C02(1:N-2,:);
A=A(1:N-2,:);
D2=D2(1:N-2,:);
coF=coF(1:N-2,:);
LHSkron1 = kron(C02,A)+kron(D2,coF);
U02big = kron(C02,C02);
Rbig = kron(C02,coF);

for d = 1:T-1
    for k = (1-N)/2:(N-1)/2
        RHS=coF*CFS(:,:,k+(N+1)/2,d)*C02.';
        RHS=RHS(1:N-2,1:N-2);
        RHS=reshape(RHS,N^2-4*N+4,1);
        RHS(N^2-4*N+5:N^2)=0;
        LHSkron = Rbig-alph*(LHSkron1-k^2*U02big)/T;
        for j=1:N
            LHSkron(N^2-4*N+4+j,:)=kron(ones(1,N),[zeros(1,(j-1)) 1 zeros(1,N-j)]);
            LHSkron(N^2-3*N+4+j,:)=kron((-1).^linspace(1,N,N),[zeros(1,(j-1)) 1 zeros(1,N-j)]);
        end
        for j=1:N-2
            LHSkron(N^2-N+2+j,:)=[zeros(1,(j-1)*N) (-1).^linspace(1,N,N) zeros(1,N^2-j*N)];
            LHSkron(N^2-2*N+4+j,:)=[zeros(1,(j-1)*N) ones(1,N) zeros(1,N^2-j*N)];
        end
        Kapprox = LHSkron\RHS;
        CFS(:,:,k+(N+1)/2,d+1) = reshape(Kapprox,N,N);%bartelsStewart(A-k^2*U02,U02.',coF,D2.',-RHS);%reshape(LHS\RHS,[],N);
    end
end

kapproxl=zeros(N,N,N,d);
for d = 1:T
    FFapprox = reshape(CFS(:,:,:,d),N,N^2);
    FFapprox = chebtech2.coeffs2vals(FFapprox);
    FFapprox = reshape(FFapprox,N,N,N);
    FFapprox = permute(FFapprox,[2 1 3]);
    FFapprox = reshape(FFapprox,N,N^2);
    FFapprox = chebtech2.coeffs2vals(FFapprox);
    FFapprox = reshape(FFapprox,N,N,N);
    FFapprox = permute(FFapprox,[3 2 1]);
    FFapprox = reshape(FFapprox,N,N^2);
    FFapprox = real(trigtech.coeffs2vals(FFapprox));
    FFapprox = reshape(FFapprox,N,N,N);
    FFapprox = permute(FFapprox,[2 1 3]);
    kapproxl(:,:,:,d)=FFapprox; %t,r,z,time
end


%% Slice plot of the results

% Remove the doubled up data in the grid
tt2 = tt(floor(N/2)+1:end,:,:);
rr2 = rr(floor(N/2)+1:end,:,:);
zz2 = zz(floor(N/2)+1:end,:,:);

% Slices in the cylinder to plot
rslice = rr2(floor(N/4)+1,1,1);
zslice = squeeze(zz2(1,1,[floor(N/4)+1 floor(N/2)+1]));
tslice = tt2(1,[1 floor(N/4)+1 floor(3*N/4)+1]+2,1)';

% Add in the theta=pi to close the cylinder in the plots
tt2 = [tt2 pi*ones(floor(N/2)+1,1,N)];
rr2 = [rr2 rr2(:,1,:)];
zz2 = [zz2 zz2(:,1,:)];

for d = 1:T
    ff2 = reshape(kapproxl(:,:,:,d),[N N N]);
    % Remove the doubled up data in the grid
    ff2 = ff2(floor(N/2)+1:end,:,:);
    ff2 = [ff2 repmat(ff2(:,1),[1 1 N])];
    hslicer = slice(tt2,rr2,zz2,ff2,tslice,rslice,zslice);
    hold on
    for j = 1:numel(hslicer)
        h = hslicer(j);
        [xs,ys,zs] = pol2cart(h.XData,h.YData,h.ZData);
        surf(xs,ys,zs,h.CData,'EdgeColor','none','FaceColor','Interp');
        hold on
    end
    delete(hslicer);
    axis([-1 1 -1 1 -1 1]);
    daspect([1 1 1])
    hold off
    % Add lighting and shine
    camlight; lighting phong
    % Change the colormap
    colormap(jet)
    title(sprintf('time = %1.3f',d*1/T))
    colorbar;
    caxis([-.05 .05]);
    pause;
end


%%
function [C1, E] = zeroDOF(C1, C2, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    for kk = 1:size(C1, 1)
        if ( abs(C1(kk,ii)) > 10*eps )
            c = C1(kk, ii); % Constant required to zero entry out.
            C1(kk,:) = C1(kk,:) - c*B(ii,:);
            E(kk,:) = E(kk,:) - c*G(ii,:)*C2.';
        end
    end
end

end
