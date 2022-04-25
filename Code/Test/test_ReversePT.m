function test_ReversePT( N, Ops )
%% Test reversePT
f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
g = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

FF = map2grid_cyl(f,N,'rzt');
GG = map2grid_cyl(g,N,'rzt');

[Fr, Ft, Fz] = reversePT(FF,GG);% Reverse without multiplication by R^2
%[Fr2, Ft2, Fz2] = reversePT2(FF,GG);
%plot_cyl(C2V_cyl(Fr,'rzt'),1,'rzt')
%Fr = mult_cyl(Ops.R0^2,Fr);
%Ft = mult_cyl(Ops.R0^2,Ft);
%Fz = mult_cyl(Ops.R0^2,Fz);
assert( infnorm_cyl(div_cyl(Fr,Fz,Ft)) < N*eps )
disp(' ''ReversePT'' works!')
end

