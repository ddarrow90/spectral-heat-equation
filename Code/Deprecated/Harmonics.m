function [HarmonicTL, HarmonicTR, HarmonicPL, HarmonicPR] = Harmonics(N)

% Test example: 
u = @(t,r,z) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^2);
%f = @(t,r,z) 2*r.^2.*(r.^2+6*z.^2-7).*sin(t).*cos(t);

% Discretization size: 
%N = 50;

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

Nold = N;
N = ceil(1.5*N);

% Construct useful ultraspherical operators: 
% Conversion operator: (ChebT -> ChebU)
U01 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); 
U01(1, 1) = 1; 

% Conversion operator: (ChebU -> C^2)
K = 1./(1:N)';
U12 = spdiags([K -K], [0 2], N, N);

% Conversion operator: (ChebT -> C^2)
U02 = U12*U01;

% First-order diff: (ChebT -> ChebU)
D1 = spdiags((0:N)', 1, N, N); 

% Second-order diff: (ChebT -> C^2)
D2 = spdiags(2*(0:N)', 2, N, N); 

% Multiplication by "r": (in C^2)
K = (1:N)'./(4:2:2*N+2)';
K1 = (3:N+2)'./(4:2:2*N+2)';
R = spdiags([K K1], [-1 1], N, N); 

% Construction of "r"-part of Laplacian on cylinder:
A = R^2*D2 + R*U12*D1;   % A[u] = r^2*u_rr + r*u_r

% Construction of "f"-part of Laplacian on cylinder
coF = R^2*U02;
% Kapprox = zeros(N^2,1);
%U02 = U02(1:N-2,:);
%A = A(1:N-2,:);
%D2 = D2(1:N-2,:);
%coF = coF(1:N-2,:);
% LHSkron1 = kron(U02, A) + kron(D2, coF);
% U02big = kron(U02, U02);

% BCs = [kron(ones(1,N), eye(N)) ;
%      kron((-1).^(1:N), eye(N)) ; 
%      kron( eye(N-2,N), ones(1,N) ) ;
%      kron( eye(N-2,N), (-1).^(1:N) ) ];

N = Nold;
A = A(1:N,1:N);
U01 = U01(1:N,1:N);
U02 = U02(1:N,1:N);
R = R(1:N,1:N);
D1 = D1(1:N,1:N);
D2 = D2(1:N,1:N);
coF = coF(1:N,1:N);
bc = [ (-1).^(0:N-1) > 0 ; 
       (-1).^(0:N-1) < 0 ]; 

%Harmonics = zeros(N,N,N,3*N^2);
for h = -(N-1)/2 : (N-1)/2
    for j = 1:3*N
        (h+(N+1)/2-1)*3*N+j
        Fapprox = zeros(N,N,N);
        if j<= N
            bc2 = [ zeros(1,j-1) 1 zeros(1,N-j); zeros(1,j-1) 0 zeros(1,N-j) ];
            bc3 = .5*[ 1-(-1)^j zeros(1,N-1); 1+(-1)^j zeros(1,N-1) ];
            bcheck = 1; 
            % r = 1: bc*X = bc2
            % z = \pm 1: X*bc.' = bc3.'
            % f(r,z,th) = T_{j-1}(z)*exp(i*h*th) 
        elseif j <= 2*N
            bc2 = [ zeros(1,j-N-1) 1 zeros(1,2*N-j); zeros(1,j-N-1) 0 zeros(1,2*N-j) ];
            bc3 = .5*[ 1+(-1)^j zeros(1,N-1); 1-(-1)^j zeros(1,N-1) ];
            bcheck = 2; 
            % z = \pm 1: X*bc.' = bc2.'
            % r = 1: bc*X = bc3
            % f(r,z,th) = T_{j+N-1}(r)*exp(i*h*th); Continuity at center?
        else
            bc2 = [ zeros(1,j-2*N-1) 0 zeros(1,3*N-j); zeros(1,j-2*N-1) 1 zeros(1,3*N-j) ];
            bc3 = .5*[ 0 1-(-1)^j zeros(1,N-2); 0 1+(-1)^j zeros(1,N-2) ];
            bcheck = 3; 
            % z = \pm 1: X*bc.' = bc2.'
            % r = 1: bc*X = bc3
            % f(r,z,th) = T_{j+2*N-1}(r)*exp(i*h*th); Continuity at center?
            
        end
        
            AA = A - h^2*U02;
            BB = U02; 
            CC = coF; 
            DD = D2; 
            RHS = zeros(N);
            %RHS = RHS(1:N-2, 1:N-2);
            
                if bcheck == 1
                    [AA, RHS] = zeroDOF(AA, BB, RHS, bc, bc2);
                    [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, bc3);
                    RHS = RHS.';
                    [CC, RHS] = zeroDOF(CC, DD, RHS, bc, bc2);
                    [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, bc3);
                    RHS = RHS.';
                elseif bcheck == 2
                    [AA, RHS] = zeroDOF(AA, BB, RHS, bc, bc3);
                    [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, bc2);
                    RHS = RHS.';
                    [CC, RHS] = zeroDOF(CC, DD, RHS, bc, bc3);
                    [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, bc2);
                    RHS = RHS.';
                else
                    [AA, RHS] = zeroDOF(AA, BB, RHS, bc, bc3);
                    [BB, RHS] = zeroDOF(BB, AA, RHS.', bc, bc2);
                    RHS = RHS.';
                    [CC, RHS] = zeroDOF(CC, DD, RHS, bc, bc3);
                    [DD, RHS] = zeroDOF(DD, CC, RHS.', bc, bc2);
                    RHS = RHS.';
                end
            
            AA = AA(1:end-2,3:end);
            BB = BB(1:end-2,3:end);
            CC = CC(1:end-2,3:end);
            DD = DD(1:end-2,3:end);
            RHS = RHS(1:end-2,1:end-2);
            
            X22 = bartelsStewart(AA, BB, CC, DD, RHS );
            
            %RHS
                if bcheck == 1
                    X12 = bc2(1:2,3:end)-bc(1:2,3:end)*X22;
                    X21 = (bc3(1:2,3:end)-bc(1:2,3:end)*X22.').';
                    X11 = bc2(1:2,1:2)-bc(1:2,3:end)*X21;
                else
                    X12 = bc3(1:2,3:end)-bc(1:2,3:end)*X22;
                    X21 = (bc2(1:2,3:end)-bc(1:2,3:end)*X22.').';
                    X11 = (bc2(1:2,1:2)-bc(1:2,3:end)*X12.').';
                end
            
            X = [ X11 X12 ; X21 X22 ]; 
            %Fapprox(:, :, h+(N+1)/2) = X;
            %Harmonics(:, :, :, (h+(N+1)/2-1)*3*N+j) = Fapprox;
            
            %{
        FFapprox = reshape(Fapprox,N,N^2);
        FFapprox = chebtech1.coeffs2vals(FFapprox);
        FFapprox = reshape(FFapprox,N,N,N);
        FFapprox = permute(FFapprox,[2 1 3]);
        FFapprox = reshape(FFapprox,N,N^2);
        FFapprox = chebtech1.coeffs2vals(FFapprox);
        FFapprox = reshape(FFapprox,N,N,N);
        FFapprox = permute(FFapprox,[3 2 1]);
        FFapprox = reshape(FFapprox,N,N^2);
        FFapprox = real(trigtech.coeffs2vals(FFapprox));
        FFapprox = reshape(FFapprox,N,N,N);
        FFapprox = permute(FFapprox,[2 1 3]);
        %}
        
        HarmonicTL = zeros(N^3,3*N^2);
        HarmonicTR = zeros(N^3,3*N^2);
        HarmonicPL = zeros(N^3,3*N^2);
        HarmonicPR = zeros(N^3,3*N^2);
        
        % Construct integration operators
        ONESval = chebtech2.coeffs2vals(ones(N,1)).';
        CC = [2 kron(2./(1-(2:2:N).^2),[0 1])];
        CC2 = kron(1./((4:4:N+3)-2),[0 1 0 -1]);
        CC2 = CC2(1:N);
        CC2 = CC2 + .5*CC;
        CCval = chebtech2.coeffs2vals(CC.').';
        CC2val = chebtech2.coeffs2vals(CC2.').';
        A = 1+(-1).^(1:N);
        Aval = chebtech2.coeffs2vals(A.').';
        DFT = trigtech.coeffs2vals(eye(N));
        DCT = chebtech2.coeffs2vals(eye(N));
        M = kron(DFT,kron(DCT,DCT));
        
        R0 = U02\R;
        D0 = U01\D1;
        Dt = spdiags(-(((1-N)/2):((N-1)/2)).^2.',0,N,N);
        
        % Construct R,T,Z components of vector harmonic, with the scalar
        % used as a Toroidal field. Z=0 here.
        Temp = zeros(N,N,N);
        Temp3 = zeros(N,N,N);
        Hr = zeros(N,N,N);
        Ht = zeros(N,N,N);
        Hr2 = zeros(N,N,N);
        Ht2 = zeros(N,N,N);
        Hz2 = zeros(N,N,N);
        Rinv = kron(1./chebpts( N ),ones(1,N));
        %{
        for jj = -(N-1)/2 : (N-1)/2
            Temp2 = D0*1i*j*Fapprox(:, :, jj+(N+1)/2);
            Temp(:,:,jj+(N+1)/2) = chebtech2.coeffs2vals(Temp2);
            Hr(:,:,jj+(N+1)/2) = Rinv.*chebtech2.coeffs2vals(1i*jj.*Fapprox(:, :, jj+(N+1)/2));
        end
        %}
        Temp = D0*1i*h*X;
        Temp = chebtech2.coeffs2vals(Temp);
        Hr(:,:,h+(N+1)/2) = Rinv.*chebtech2.coeffs2vals(chebtech2.coeffs2vals(1i*h*X).').';
        Hr((N+1)/2,:,h+(N+1)/2)=Temp((N+1)/2,:);
        %{
        for jj = 1:N
            Hr(:,:,jj) = chebtech2.vals2coeffs(Hr(:,:,jj));
        end
        for jj = 1:N
            Ht(:,:,jj) = -U01\D1*Fapprox(:, :, jj);
        end
        %}
        Hr(:,:,h+(N+1)/2) = chebtech2.vals2coeffs(chebtech2.vals2coeffs(Hr(:,:,h+(N+1)/2)).').';
        Ht(:,:,h+(N+1)/2) = -U01\D1*X;
        
        
        % Construct R,T,Z components of vector harmonic, with the scalar
        % used as a Poloidal field.
        
        
        for jj = 1:N
            Hr2(:,:,jj) = -Ht(:, :, jj)*D0.';
            Ht2(:,:,jj) = Hr(:, :, jj)*D0.';
        end
        
        
        for jj = -(N-1)/2 : (N-1)/2
            Hz2(:,:,jj+(N+1)/2) = D0*Ht(:, :, jj+(N+1)/2);
            Hz2(:,:,jj+(N+1)/2) = Hz2(:,:,jj+(N+1)/2)-1i*jj.*Hr(:,:,jj+(N+1)/2);
        end
        for jj = 1:N
            Temp2 = D0*Ht(:, :, jj);
            Temp(:,:,jj) = chebtech2.coeffs2vals(Temp2);
            Temp3(:,:,jj) = Rinv.*chebtech2.coeffs2vals(Ht(:, :, jj));
        end
        %Temp3
        Temp3((N+1)/2,:,:)=Temp((N+1)/2,:,:);
        %Temp3
        for jj = 1:N
            Temp3(:,:,jj) = chebtech2.vals2coeffs(Temp3(:,:,jj));
        end
        %Temp3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM?
        Hz2=Hz2+Temp3;
        
        % Build the things to diagonalize, for the toroidal case.
        Temp = reshape(Ht,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag1 = permute(Temp,[2 3 1]); % R, Z, Theta
        
        for jj = 1:N
            Temp(:,:,jj) = (speye(N)+R0*D0)*Ht(:,:,jj);
        end
        for jj = -(N-1)/2 : (N-1)/2
            Temp3(:,:,jj+(N+1)/2) = -1i*jj.*Hr(:, :, jj+(N+1)/2);
        end
        Temp = Temp + Temp3;
        
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag2 = permute(Temp,[2 3 1]);
        
        for jj = 1:N
            %Temp(:,:,jj) = U01\(D1*Hz(:,:,jj));
            Temp3(:,:,jj) = -Hr(:,:,jj)*D0.';
        end
        Temp = Temp3; %Temp + Temp3;
        
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag3 = permute(Temp,[2 3 1]);
        
        %{
        Temp = reshape(Hz,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag4 = permute(Temp,[2 3 1]);
        %}
        diag4=zeros(N,N,N);
        
        %diag5 = diag1;
        
        for jj = 1:N
            Temp(:,:,jj) = R0*Hr(:,:,jj);
        end
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag6 = permute(Temp,[2 3 1]);
        
        % Kronecker-ing everything up.
        diag1 = reshape(diag1,N^3,1).';
        diag2 = reshape(diag2,N^3,1).';
        diag3 = reshape(diag3,N^3,1).';
        diag4 = reshape(diag4,N^3,1).';
        diag5 = diag1;
        diag6 = reshape(diag6,N^3,1).';
        
        % Build Kronecker-ed boundary rows for the toroidal case.
        B1 = kron(ones(1,N),kron(CCval,ONESval)).*diag1;
        B1 = B1 + kron(ones(1,N),kron(CCval,CC2val)).*diag2;
        
        HarmonicTL(:,(h+(N+1)/2-1)*3*N+j) = (B1*M).';
        
        Temp4 = kron(ones(1,N),kron(CCval,ONESval)).*diag3;
        B2 = Temp4*M;
        Temp4 = -kron(ones(1,N),kron(CCval,ONESval)).*diag4;
        B2 = B2 + Temp4*kron(DFT,kron(DCT,chebtech2.coeffs2vals(D0*eye(N))));
        Temp4 = kron(ones(1,N),kron(Aval,CC2val)).*diag5;
        B2 = B2 + Temp4*kron(trigtech.coeffs2vals(Dt*eye(N)),kron(DCT,DCT));
        Temp4 = kron(ones(1,N),kron(Aval,CC2val)).*diag6;
        B2 = B2 + Temp4*kron(DFT,kron(DCT,chebtech2.coeffs2vals(D0*eye(N))));
        HarmonicTR(:,(h+(N+1)/2-1)*3*N+j) = B2.';

        % Build the things to diagonalize, for the poloidal case.
        
        Temp = reshape(Ht2,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag1 = permute(Temp,[2 3 1]); % R, Z, Theta
        
        for jj = 1:N
            Temp(:,:,jj) = (speye(N)+R0*D0)*Ht2(:,:,jj);
        end
        for jj = -(N-1)/2 : (N-1)/2
            Temp3(:,:,jj+(N+1)/2) = -1i*jj.*Hr2(:, :, jj+(N+1)/2);
        end
        Temp = Temp + Temp3;
        
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag2 = permute(Temp,[2 3 1]);
        
        for jj = 1:N
            Temp(:,:,jj) = D0*Hz2(:,:,jj);
            Temp3(:,:,jj) = -Hr2(:,:,jj)*D0.';
        end
        Temp = Temp + Temp3;
        
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag3 = permute(Temp,[2 3 1]);
        
        Temp = reshape(Hz2,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag4 = permute(Temp,[2 3 1]);
        
        %diag5 = diag1;
        
        for jj = 1:N
            Temp(:,:,jj) = R0*Hr2(:,:,jj);
        end
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[2 1 3]);
        Temp = reshape(Temp,N,N^2);
        Temp = chebtech1.coeffs2vals(Temp);
        Temp = reshape(Temp,N,N,N);
        Temp = permute(Temp,[3 2 1]);
        Temp = reshape(Temp,N,N^2);
        Temp = real(trigtech.coeffs2vals(Temp));
        Temp = reshape(Temp,N,N,N);
        diag6 = permute(Temp,[2 3 1]);
        
        % Kronecker-ing everything up.
        diag1 = reshape(diag1,N^3,1).';
        diag2 = reshape(diag2,N^3,1).';
        diag3 = reshape(diag3,N^3,1).';
        diag4 = reshape(diag4,N^3,1).';
        diag5 = diag1;
        diag6 = reshape(diag6,N^3,1).';
        
        % Build Kronecker-ed boundary rows for the poloidal case.
        B1 = kron(ones(1,N),kron(CCval,ONESval)).*diag1;
        B1 = B1 + kron(ones(1,N),kron(CCval,CC2val)).*diag2;
        
        HarmonicPL(:,(h+(N+1)/2-1)*3*N+j) = (B1*M).';
        
        Temp4 = kron(ones(1,N),kron(CCval,ONESval)).*diag3;
        B2 = Temp4*M;
        Temp4 = -kron(ones(1,N),kron(CCval,ONESval)).*diag4;
        B2 = B2 + Temp4*kron(DFT,kron(DCT,chebtech2.coeffs2vals(D0*eye(N))));
        Temp4 = kron(ones(1,N),kron(Aval,CC2val)).*diag5;
        B2 = B2 + Temp4*kron(trigtech.coeffs2vals(Dt*eye(N)),kron(DCT,DCT));
        Temp4 = kron(ones(1,N),kron(Aval,CC2val)).*diag6;
        B2 = B2 + Temp4*kron(DFT,kron(DCT,chebtech2.coeffs2vals(D0*eye(N))));
        HarmonicPR(:,(h+(N+1)/2-1)*3*N+j) = B2.';
    end  
end

HarmonicTL=1/N^3*sparse(HarmonicTL);
HarmonicTR=1/N^3*sparse(HarmonicTR);
HarmonicPL=1/N^3*sparse(HarmonicPL);
HarmonicPR=1/N^3*sparse(HarmonicPR);
end

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