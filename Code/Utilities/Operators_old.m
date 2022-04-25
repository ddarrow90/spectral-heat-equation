classdef Operators
    properties
        C01 % Convert Cheb-T to Cheb-U
        C12 % Convert Cheb-U to C^2
        C02 % Convert Cheb-T to C^2
        D0 % Differentiate Cheb-T once
        D1 % Differentiate Cheb-T once and convert to Cheb-U
        D2 % Differentiate Cheb-T twice and convert to C^2
        R0 % Multiply Cheb-T by active variable
        R2 % Multiply C^2 by active variable
        A % R2^2*D2 + R2*C12*D1. Used in Laplacian.
        coF % R2^2*C02. Used in Laplacian.
        Dt % Differentiate Fourier once
        DCT % Values at Chebyshev nodes to Cheb-T coefficients
        DFT % Values at equidistant points to Fourier coefficients
        iDCT % Cheb-T to values at Chebyshev nodes
        iDFT % Fourier to values at equidistant points
        CC % Integrate Cheb-T over [-1, 1]
        CC2 % Integrate Cheb-T over [0, 1]
        CCf % Integrate Fourier over [0, 2*pi]
    end
    methods
        function obj = Operators(N,arg)
            if ( ~isnumeric(N) )
                error('Must be an integer value!')
            elseif ( floor(N) ~= N )
                error('Must be an integer value!')
            else
                obj.C01 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); 
                obj.C01(1, 1) = 1; 
                
                K = 1./(1:N)';
                obj.C12 = spdiags([K -K], [0 2], N, N);
                
                obj.C02 = obj.C01*obj.C12;
                
                obj.D1 = spdiags((0:N)', 1, N, N);
                
                obj.D2 = spdiags(2*(0:N)', 2, N, N);
                
                obj.D0 = obj.C01\obj.D1;
                
                obj.R0 = spdiags([1 .5; .5*ones(N,1) .5*ones(N,1)],[-1,1],N,N);
                
                K = (1:N)'./(4:2:2*N+2)';
                K1 = (3:N+2)'./(4:2:2*N+2)';
                obj.R2 = spdiags([K K1], [-1 1], N, N);
                
                obj.A = obj.R2^2*obj.D2 + obj.R2*obj.C12*obj.D1;
                
                obj.coF = obj.R2^2*obj.C02;
                
                if mod(N,2)
                	obj.Dt = spdiags((1i*((1-N)/2:(N-1)/2))', 0, N, N); 
                else
                	obj.Dt = spdiags((1i*((-N)/2:N/2-1))', 0, N, N); 
                end
                
                obj.iDCT = chebtech2.coeffs2vals(eye(N));
                
                obj.iDFT = trigtech.coeffs2vals(eye(N)); 
                
                obj.DCT = chebtech2.vals2coeffs(eye(N)); 
                
                obj.DFT = trigtech.vals2coeffs(eye(N)); 
                
                obj.CC = [2 kron(2./(1-(2:2:N).^2),[0 1])];
                obj.CC = obj.CC(1:N);
                
                obj.CC2 = kron(1./((4:4:N+3)-2),[0 1 0 -1]); % Integrate on [0,1]
                obj.CC2 = obj.CC2(1:N);
                obj.CC2 = obj.CC2 + .5*obj.CC;
                obj.CCf = [zeros(1,(ceil((N-1)/2))) 2*pi zeros(1,(floor((N-1)/2)))];
            end
        end
    end
end