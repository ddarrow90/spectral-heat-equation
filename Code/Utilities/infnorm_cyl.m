function m = infnorm_cyl(varargin)
% This function finds the maximum value within a collection of arrays, of 
% any size and any number of dimensions.

% Dave Darrow. February 3, 2018.

v = zeros(nargin,1);
for i = 1:nargin
    v(i) = norm(varargin{i}(:),inf);
end
m = norm(v,inf);

end