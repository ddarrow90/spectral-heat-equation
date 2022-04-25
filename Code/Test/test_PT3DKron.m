function test_PT3DKron( N, Ops )
%% Test PT3DKron
f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
g = @(r,z,t) 0*r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

FF = map2grid_cyl(f,N);
GG = map2grid_cyl(g,N);

[Fr, Ft, Fz] = reversePT(FF,GG);% Reverse without multiplication by R

DIV = div_cyl(Fr,Ft,Fz);
assert( norm(DIV(:),inf) < N*eps )
disp(' ''reversePT'' works!')
[T, P] = PT3DKron(Fr,Ft,Fz,N);

[Fr2, Ft2, Fz2] = reversePT(T,P);

temp = infnorm_cyl(Fr-Fr2,Ft-Ft2,Fz-Fz2);
assert(temp < N*eps)
disp(' ''PT3DKron'' works!')
% The reason this is off is because of expected errors in the reversePT
% code. Because we do not switch bases (and then divide by R), we incur an
% error in the highest frequency mode.

end

