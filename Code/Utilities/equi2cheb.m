function vals = equi2cheb(vals,a)
% This code takes as argument a matrix of function values at N-dimensional
% Cartesian products of finite sets of points, where some of the N
% dimensions are equally sampled. It resamples these dimensions (according
% to the second argument, a string of boolean values) at Chebyshev points.
% 
% In the "order" string, the characters 't', 'T', or '1' denote 
% "trigonometric" and equispaced dimensions respectively. These are the 
% ones for which the DFT is used. Any other character is used to denote a 
% dimension not to transform.
% 
% September 17, 2017. David Darrow.

if ( size(a,1) > 1 )
    return
end

dim = size(a,2);
SIZE = size(vals);

if isstring(a)
    a = char(a);
end

for kk = 1:dim
    vals = reshape(vals,SIZE(kk),prod(SIZE)/SIZE(kk));
    if ( a(kk) == 't' || a(kk) == 'T')
        for ii = 1:prod(SIZE)/SIZE(kk)
            temp = chebfun(vals(:,ii),'trig');
            vals(:,ii) = sample(temp,SIZE(kk));
        end
    elseif ( a(kk) == '1')
        for ii = 1:prod(SIZE)/SIZE(kk)
            temp = chebfun(vals(:,ii),'equi');
            vals(:,ii) = sample(temp,SIZE(kk));
        end
    end

    vals = reshape(vals,circshift(SIZE,dim-kk+1));
    vals = permute(vals,circshift(1:dim,dim-1));
end

end