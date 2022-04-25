function v = LHopital_cheb(v)
% This function takes an array, wherein the first dimension represents
% coefficients of a Chebyshev T series. It attempts to divide the array by
% the variable in the first dimension, replacing the zero value by the
% derivative in the first dimension.
%
% August 22, 2017. David Darrow.

SPcheck = 0;
if issparse(v)
    v = full(v);
    SPcheck = 1;
end

SIZE = size(v);
N = SIZE(1);

if mod(N+1,2) % You don't need this function if N is even.
    return
end

if ( size(SIZE,2) > 2 )
    v = reshape(v,N,prod(SIZE)/N); % Get into matrix form.
end
 
D0 = spdiags(.5*ones(N,1)*[1 -1], [0 2], N, N); % Conversion to Cheb-U.
D0(1, 1) = 1; % Still conversion to Cheb-U.
D0 = D0\spdiags((0:N)', 1, N, N); % Differentiation in Cheb T.

Rinv = 1./chebpts(N);

for kk = 1:size(v,2)
    temp = chebtech2.coeffs2vals(D0*v(:,kk));
    temp = temp((N+1)/2);
    v(:,kk) = Rinv.*chebtech2.coeffs2vals(v(:,kk));
    v((N+1)/2,kk) = temp;
    % Note: 1/x is interpolated when appropriate.
    v(:,kk) = chebtech2.vals2coeffs(v(:,kk));
end

if ( size(SIZE,2) > 2 )
    v = reshape(v,SIZE); % Get into original form.
end

if SPcheck
    v = sparse(v);
end
end