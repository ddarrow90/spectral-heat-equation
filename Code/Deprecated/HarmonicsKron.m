function harmonics = HarmonicsKron(N)

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
coF=R^2*U02;
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
U02 = U02(1:N,1:N);
D2 = D2(1:N,1:N);
coF = coF(1:N,1:N);
% bc = [ (-1).^(0:N-1) > 0 ; 
%        (-1).^(0:N-1) < 0 ]; 
   
jj = 1; 
for h = -(N-1)/2 : (N-1)/2
%     for j = 1% 1:2:N
%         Fapprox = zeros(N,N,N);
%         if ( j <= N )
            % r = 1: bc*X = bc2
            % z = \pm 1: X*bc.' = bc3.'
            % f(r,z,th) = T_{j-1}(z)*exp(i*h*th) 
            
%             bc2 = [ zeros(1,j-1) 1 zeros(1,N-j); 
%                     zeros(1,j-1) 0 zeros(1,N-j) ];
%             bc2 = [ 1 zeros(1,N-1) ; 
%                      zeros(1, N)   ]; 
%             bc2 = [ zeros(1,j-1) -1 0 1 zeros(1,N-j-2); 
%                     zeros(1,j-1) 0 zeros(1,N-j) ];
                
%             bc3 = .5*[1-(-1)^j      zeros(1,N-1); 
%                       1+(-1)^j  zeros(1,N-1) ];
%             bc3 = [ zeros(1,j-1) 1 zeros(1,N-j); 
%                     zeros(1,j-1) 0 zeros(1,N-j) ];
%             bc3 = .5*[   0      zeros(1,N-1); 
%                          0      zeros(1,N-1) ];     
                  
%             bcheck = 1; 
            
%         elseif j <= 2*N
%             bc2 = [ zeros(1,j-N-1) 1 zeros(1,2*N-j); zeros(1,j-N-1) 0 zeros(1,2*N-j) ];
%             bc3 = .5*[ 1+(-1)^j zeros(1,N-1); 1-(-1)^j zeros(1,N-1) ];
%             bcheck = 2; 
            % z = \pm 1: X*bc.' = bc2.'
            % r = 1: bc*X = bc3
            % f(r,z,th) = T_{j+N-1}(r)*exp(i*h*th); Continuity at center?
%         else
%             bc2 = [ zeros(1,j-2*N-1) 0 zeros(1,3*N-j); zeros(1,j-2*N-1) 1 zeros(1,3*N-j) ];
%             bc3 = .5*[ 0 1-(-1)^j zeros(1,N-2); 0 1+(-1)^j zeros(1,N-2) ];
%             bcheck = 3; 
            % z = \pm 1: X*bc.' = bc2.'
            % r = 1: bc*X = bc3
            % f(r,z,th) = T_{j+2*N-1}(r)*exp(i*h*th); Continuity at center?          
%         end
        
            AA = A - h^2*U02;
            BB = U02; 
            CC = coF; 
            DD = D2; 
            
%             AA = AA(1:end-2,:);
%             BB = BB(1:end-2,:);
%             CC = CC(1:end-2,:);
%             DD = DD(1:end-2,:); 
    
            [U, S, V] = svd(full(kron(BB, AA) + kron(DD, CC)));
%             idx = find(diag(S)./S(1,1)>1e-7,1,'last');
%             [N^2 - idx + 1, 2*N-4]
            for j = (N^2-2*N+4)+1:N^2
                X = reshape(V(:,j), N, N); 
            
%             BC1 = kron(speye(N),bc);
%             BC2 = kron(bc,speye(N));
                  
%             bc3 = bc3.'; 
%             X = [ kron(BB, AA) + kron(DD, CC) ; 
%                         BC1(1:end-2,:)        ;
%                         BC2(1:end-2,:)        ] \ ...
%                                           [ zeros((N-2)^2,1)  ; 
%                                               bc2(1:end-2).'  ; 
%                                               bc3(1:end-2).' ];             
%             X = reshape(X, N, N);
%             G = chebtech2.coeffs2vals(chebtech2.coeffs2vals(X).').';
%             g = chebfun2( G ); 
%             plot(g)
            
%             gg = chebfun2(@(z,r) besselj(h, j*r).*exp(-j*z) );
%             X = coeffs2( gg, N, N); 
%             Fapprox(:, :, h+(N+1)/2) = X;
                harmonics{jj} = X;
                jj = jj + 1; 
            
                % residue check:
            AA = A - h^2*U02;
            BB = U02; 
            CC = coF; 
            DD = D2; 
            norm(AA*X*BB.'+CC*X*DD.')
            end
end

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