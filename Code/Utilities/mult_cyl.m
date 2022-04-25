function PROD = mult_cyl(A,X,B)
% This code takes two matrices, A and B, as well as a 3D array X (with the 
% dimension of each X(:,:,k) such that A*X(:,:,k)*B.' is defined). It
% outputs the array PROD with PROD(:,:,k) = A*X(:,:,k)*B.'. If only two
% arguments are given, it performs the same calculation but without A, B, 
% or X, depending on the nature of the arguments.

% Dave Darrow. February 3, 2018.



if ( nargin == 3 ) % Matrices on either side of a 3D array.
    assert( ismatrix(A) && ndims(X) == 3 && ismatrix(B) );
    
    PROD = zeros(size(A,1),size(B,1),size(X,3));
    for i = 1:size(X,3)
        PROD(:,:,i) = A*X(:,:,i)*B.';
    end
else
    if ( ismatrix(A) && ndims(X) == 3 ) % Matrix on left of 3D array.
        PROD = zeros(size(A,1),size(X,2),size(X,3));
        for i = 1:size(X,3)
            PROD(:,:,i) = A*X(:,:,i);
        end
    elseif ( ndims(A) == 3 && ismatrix(X) ) % Matrix on right of 3D array.
        PROD = zeros(size(A,1),size(X,1),size(A,3));
        for i = 1:size(A,3)
            PROD(:,:,i) = A(:,:,i)*X.';
        end
    elseif ( ismatrix(A) && ismatrix(X) ) % Two matrices.
        PROD = A*X.';
    else
        error('The function mult_cyl requires different arguments.')
    end
end


end