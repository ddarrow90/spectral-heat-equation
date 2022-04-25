function [H, Indices] = RemoveDependencies(H,tol)
% This code takes as input some matrix and some tolerance; by taking the QR
% decomposition of this matrix and sorting the eigenvalues of the "R"
% matrix by size, it determines the linearly dependent columns of the input
% matrix. It then outputs a matrix with only the linearly independent
% columns. The second output is a column vector with the indices of the 
% preserved rows of the input matrix.
%
% August 25, 2017. Alex Townsend, David Darrow.

[~, R, E] = qr(H.');
% Harmonic.'*E = Q*R.
if ( nargin < 2 )
    tol = 1e-14;
end
idx = find(abs(diag(R))>tol,1,'last');
Indices = find(all(E(:,idx+1:size(H,1)).' == 0)).';

% Remove the linearly dependent rows.
H = H(Indices,:);

end