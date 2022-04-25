function row = L2IPR(Hr,Hz,Ht)
% L^2 Inner Product Row.
% 
% This code takes as input the R, Z, and Theta components of some
% cylindrical function H, in the form of coefficient arrays for
% Chebyshev-Chebyshev-Fourier series. The order of variables (in terms of
% dimensions of the coefficient arrays) is assumed to be R, Z, and then 
% Theta. 
%
% The program produces a row vector that, when right multiplied by an 
% argument function F, gives the L^2 inner product (in the cylinder) 
% between H and F. F must be in the form [T; P], where T and P are the
% (vectorized) toroidal and poloidal components of F (with the same order
% of variables as before).
%
% The cylinder is assumed to be [0, 1]x[-1, 1]x[0, 2*pi] for R, Z, and
% Theta respectively.
%
% August 23, 2017. David Darrow.

if ( ~all(size(Hr)==size(Hz)) || ~all(size(Hr)==size(Ht)) )
    error('Discretization sizes must be the same between vector components')
end

N1 = size(Hr,1); N2 = size(Hr,2); N3 = size(Hr,3);
Ntot = N1*N2*N3;

opsR = Operators(N1,'D0 iDCT DCT R0 CC2');
Dr = opsR.D0;
iDCTr = opsR.iDCT;
DCTr = opsR.DCT;
R = opsR.R0;
CCr = opsR.CC2;

opsZ = Operators(N2,'D0 iDCT DCT CC');
Dz = opsZ.D0;
iDCTz = opsZ.iDCT;
DCTz = opsZ.DCT;
CCz = opsZ.CC;

opsT = Operators(N3,'Dt iDFT DFT CCf');
Dt = opsT.Dt;
iDFT = opsT.iDFT;
DFT = opsT.DFT;
CCt = opsT.CCf;

M = kron(iDFT,kron(iDCTz,iDCTr));
DrBIG = kron(speye(N3),kron(speye(N2),Dr));
DtBIG = kron(Dt,kron(speye(N2),speye(N2)));

Alt = 1+(-1).^(1:N2);

% Build the things to diagonalize, for the toroidal case.
diag1 = C2V_cyl(Ht, 'rzt');
diag1 = reshape(diag1,1,Ntot);

diag2 = zeros(N1,N2,N3);
for h = (1-N3)/2:(N3-1)/2 % Ht + R*Ht_r - Ht_theta
    diag2(:,:,h+(N3+1)/2) = (speye(N1)+R*Dr)*Ht(:,:,h+(N3+1)/2)-1i*h.*Hr(:, :, h+(N3+1)/2);
end
diag2 = C2V_cyl(diag2, 'rzt');
diag2 = reshape(diag2,1,Ntot);

diag3 = zeros(N1,N2,N3);
for h = (1-N3)/2:(N3-1)/2 % Hr_z - Hz_r
    diag3(:,:,h+(N3+1)/2) = Hr(:,:,h+(N3+1)/2)*Dz.'-Dr*Hz(:, :, h+(N3+1)/2);
end
diag3 = C2V_cyl(diag3, 'rzt');
diag3 = reshape(diag3,1,Ntot);

diag4 = C2V_cyl(Hz, 'rzt');
diag4 = reshape(diag4,1,Ntot);

diag5 = diag1;

diag6 = zeros(N1,N2,N3);
for h = (1-N3)/2:(N3-1)/2 % Hr_z - Hz_r
    diag6(:,:,h+(N3+1)/2) = R*Hr(:,:,h+(N3+1)/2);
end
diag6 = C2V_cyl(diag6, 'rzt');
diag6 = reshape(diag6,1,Ntot);

diag7 = zeros(N1,N2,N3);
for h = (1-N3)/2:(N3-1)/2 % Hr_z - Hz_r
    diag7(:,:,h+(N3+1)/2) = (speye(N1)+R*Dr)*Hr(:,:,h+(N3+1)/2)*Dz.'-(Dr+R*Dr^2)*Hz(:, :, h+(N3+1)/2);
    diag7(:,:,h+(N3+1)/2) = diag7(:,:,h+(N3+1)/2) + h^2*LHopital_cheb(Hz(:, :, h+(N3+1)/2)) + 1i*h*Ht(:, :, h+(N3+1)/2)*Dz.';
end
diag7 = C2V_cyl(diag7, 'rzt');
diag7 = reshape(diag7,1,Ntot);


% Build Kronecker-ed boundary rows for the toroidal case.
B1 = kron(CCt*DFT,kron(CCz*DCTz,ones(1,N1)*DCTr)).*diag1;
B1 = B1 + kron(CCt*DFT,kron(CCz*DCTz,CCr*DCTr)).*diag2;

B1 = B1*M;

B2 = (kron(CCt*DFT,kron(CCz*DCTz,ones(1,N1)*DCTr)).*diag3)*M;
B2 = B2 + (kron(CCt*DFT,kron(CCz*DCTz,ones(1,N1)*DCTr)).*diag4)*M*DrBIG;
B2 = B2 - (kron(CCt*DFT,kron(Alt*DCTz,CCr*DCTr)).*diag5)*M*DtBIG;
B2 = B2 - (kron(CCt*DFT,kron(Alt*DCTz,CCr*DCTr)).*diag6)*M*DrBIG;
B2 = B2 + (kron(CCt*DFT,kron(CCz*DCTz,CCr*DCTr)).*diag7)*M;

row = [B1 B2];
end