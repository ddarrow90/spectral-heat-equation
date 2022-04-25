function test_mult_cyl( N, Ops )
%% Test mult_cyl

assert( norm(mult_cyl(Ops.R0,Ops.Dt)-Ops.R0*Ops.Dt.',inf) < eps );

f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);
FF = func2grid(f,N,'rzt',1);

multFF = zeros(N,N,N);
for i = 1:N
    multFF(:,:,i) = Ops.R0*FF(:,:,i);
end
A = mult_cyl(Ops.R0,FF)-multFF;
assert( norm(A(:),inf) < eps );

multFF = zeros(N,N,N);
for i = 1:N
    multFF(:,:,i) = FF(:,:,i)*Ops.R0.';
end
A = mult_cyl(FF,Ops.R0)-multFF;
assert( norm(A(:),inf) < eps );

multFF = zeros(N,N,N);
for i = 1:N
    multFF(:,:,i) = Ops.D0*FF(:,:,i)*Ops.R0.';
end

A = mult_cyl(Ops.D0,FF,Ops.R0)-multFF;
assert( norm(A(:),inf) < eps );
disp(' ''mult_cyl'' works!')
end

