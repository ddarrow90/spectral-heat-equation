function test_LHopital_cheb( N, ~ )
%% Test LHopital_cheb
f = @(r,z,t) r.^2.*(z-1).*(1-r.*sin(t));
FF = func2grid(f,N,'rzt');
FF = V2C_cyl(FF,'rzt');
FF = LHopital_cheb(FF);
g = @(r,z,t) r.*(z-1).*(1-r.*sin(t));
GG = func2grid(g,N,'rzt');
GG = V2C_cyl(GG,'rzt');
assert( norm(FF(:)-GG(:)) < N*eps )

disp(' ''LHopital_cheb'' works!')
end

