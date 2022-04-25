function [Harmonic, NORM] = HarmonicToRows(N,tol)
% This code only needs the intended discretization size as an input. It
% first generates a basis of (discrete) cylindrical harmonics by finding
% the nullspace of the Laplacian (using the SVD). Then it proceeds to turn
% each harmonic into a row vector that calculates the harmonic's L^2 inner
% product (in the cylinder) with another function. The program outputs a
% matrix of these row vectors for use as boundary conditions.
% 
% As its final calculation, the program removes all linearly dependent rows
% from this matrix and normalizes it. The norm is available as a second
% output of this code.
%
% August 24, 2017. David Darrow.

Q = 2*N-3; % Number of harmonics for each Fourier mode.
% Two additional harmonics are added in the constant Fourier case.

% Make sure the parity of N is odd:
if ( mod(N+1, 2) )  
    N = N + 1;
end

% Create a whole bunch of operators. These will be used to construct the 
% Laplacian and to take the curl of our Toroidal and Poloidal component 
% harmonics such that we get the R, Z, and Theta components of the 
% corresponding vector.
ops = Operators(N,'C02 D0 D2 A coF');
D0 = ops.D0;
A = ops.A;
BB = ops.C02;
CC = ops.coF;
DD = ops.D2;

Harmonic = zeros(2*N*Q+4, 2*N^3); % Initialize output.
    
for h = -(N-1)/2 : (N-1)/2
    AA = A - h^2*BB;
    
	[~, ~, V] = svd(full(kron(BB, AA) + kron(DD, CC))); % The Laplacian.
%	idx = find(diag(S)./S(1,1)>1e-7,1,'last');
%	[N^2 - idx + 1, 2*N-4]
    foo = 0;
    if ( h == 0 )
        foo = 2;
    end
	for j = (N^2-Q+1-foo):N^2
        X = reshape(V(:,j), N, N);
        
        % Construct R,T,Z components of vector harmonic, with the scalar
        % used as a Toroidal field.
        
        % First, initialize everything for both Toroidal and Poloidal cases.
        Hr = zeros(N,N,N);
        Ht = zeros(N,N,N);
        Hr2 = zeros(N,N,N);
        Ht2 = zeros(N,N,N);
        Hz2 = zeros(N,N,N);
        
        % Hr = "r" part of Toroidal harmonic.
        Hr(:,:,h+(N+1)/2) = LHopital_cheb(1i*h*X);
        
        % Ht = "theta" part of Toroidal harmonic.
        Ht(:,:,h+(N+1)/2) = -D0*X;
        
        % Hz is actually zero; there is no Z component of Toroidal fields.
        
        % Construct R,T,Z components of vector harmonic, with the scalar
        % used as a Poloidal field.
        
        % Hr2 = "r" part of Poloidal harmonic.
        Hr2(:,:,h+(N+1)/2) = -Ht(:, :, h+(N+1)/2)*D0.';
        
        % Ht2 = "theta" part of Poloidal harmonic.
        Ht2(:,:,h+(N+1)/2) = Hr(:, :, h+(N+1)/2)*D0.';
        
        % Hz2 = "z" part of Poloidal harmonic.
        Hz2(:,:,h+(N+1)/2) = D0*Ht(:, :, h+(N+1)/2) + LHopital_cheb(Ht(:, :, h+(N+1)/2));
        Hz2(:,:,h+(N+1)/2) = Hz2(:,:,h+(N+1)/2) - LHopital_cheb(1i*h.*Hr(:,:,h+(N+1)/2));
        
        % Create boundary row for Toroidal case.
        row = L2IPR(Hr,zeros(N,N,N),Ht); % Again, the Z component vanishes
        if ( j == N^2-Q-1 ) % Extra cases for constant Fourier mode
            Harmonic(N*Q+1, :) = row; 
        elseif ( j == N^2-Q )
            Harmonic(N*Q+2, :) = row;
        else
            Harmonic((h+(N+1)/2-1)*Q + j+Q-N^2, :) = row;
        end

        % Create boundary row for Poloidal case.
        row = L2IPR(Hr2,Hz2,Ht2);
        if ( j == N^2-Q-1 ) % Extra cases for constant Fourier mode
            Harmonic(2*N*Q+3, :) = row;
        elseif ( j == N^2-Q )
            Harmonic(2*N*Q+4, :) = row;
        else
            Harmonic(N*Q + 2 + (h+(N+1)/2-1)*Q + j+Q-N^2, :) = row;
        end
	end
end

% Normalize the harmonic matrix. It's generally very badly scaled.
NORM = norm(Harmonic);
Harmonic = Harmonic/norm(Harmonic);
end