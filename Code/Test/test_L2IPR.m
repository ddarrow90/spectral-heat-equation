function test_L2IPR( N, Ops )
%% Test L2IPR
return;
f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
g = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

FF = func2grid(f,N,'rzt');
GG = func2grid(g,N,'rzt');
FF = V2C_cyl(FF,'rzt');
GG = V2C_cyl(GG,'rzt');
[Fr, Ft, Fz] = reversePT2(FF,GG);% Reverse without multiplication by R

Fr2 = C2V_cyl(Fr,'rzt');
Ft2 = C2V_cyl(Ft,'rzt');
Fz2 = C2V_cyl(Fz,'rzt');
[T, P, s] = PT3DKron(Fr2,Ft2,Fz2,N,1);

DotP = Fr2.^2 + Ft2.^2 + Fz2.^2;
DotP = V2C_cyl(DotP,'rzt');
magnitude = Ops.CC2*DotP(:,:,(N+1)/2)*Ops.CC.';
row = L2IPR(Fr,Fz,Ft);
row*[T(:); P(:)]
end

