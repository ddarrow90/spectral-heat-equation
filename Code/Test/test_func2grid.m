function test_func2grid( N, ~ )
%% Test func2grid

f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

FF = func2grid(f,N,'rzt');

t = pi*trigpts( N );
r = chebpts( N );

for i = 1:N
    for j = 1:N
        for k = 1:N
            assert(FF(i,j,k) - f(r(i),r(j),t(k)) < eps);
        end
    end
end

disp(' ''func2grid'' works!')
end

