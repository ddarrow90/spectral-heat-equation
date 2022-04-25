f = @(r,z,t) r.^2.*sin(t).*cos(t).*(1-r.^2).*(1-z.^4);

GroundTruth = Heat3DMatEqn(f,140,20,.01);
GroundTruth = equi2cheb(GroundTruth(:,:,:,20),'00t');
GroundTruth = chebfun3(GroundTruth);

%% New method
tic
TEST = Heat3DMatEqn_ExactSln(f,60,20,.01);
TIME(1,1) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(1,1) = norm(GroundTruth - TEST);
disp('small')

% Higher Size.
tic
TEST = Heat3DMatEqn_ExactSln(f,80,20,.01);
TIME(2,1) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(2,1) = norm(GroundTruth - TEST);
disp('medium')

% Higher size.
tic
TEST = Heat3DMatEqn_ExactSln(f,100,20,.01);
TIME(3,1) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(3,1) = norm(GroundTruth - TEST);
disp('large')

% Higher size.
tic
TEST = Heat3DMatEqn_ExactSln(f,120,20,.01);
TIME(4,1) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(4,1) = norm(GroundTruth - TEST);
disp('very large')
%% Collocation
tic
TEST = Heat3DMatEqnColloc(f,60,20,.01);
TIME(1,2) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(1,2) = norm(GroundTruth - TEST);
disp('small')

% Higher Size.
tic
TEST = Heat3DMatEqnColloc(f,80,20,.01);
TIME(2,2) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(2,2) = norm(GroundTruth - TEST);
disp('medium')

% Higher size.
tic
TEST = Heat3DMatEqnColloc(f,100,20,.01);
TIME(3,2) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(3,2) = norm(GroundTruth - TEST);
disp('large')

% Higher size.
tic
TEST = Heat3DMatEqnColloc(f,120,20,.01);
TIME(4,2) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'00t');
TEST = chebfun3(TEST);
ERROR(4,2) = norm(GroundTruth - TEST);
disp('very large')
%% Finite Difference
tic
TEST = Heat3DMatEqnFD(f,60,20,.01);
TIME(1,3) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'11t');
TEST = chebfun3(TEST);
ERROR(1,3) = norm(GroundTruth - TEST);
disp('small')

% Higher Size.
tic
TEST = Heat3DMatEqnFD(f,80,20,.01);
TIME(2,3) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'11t');
TEST = chebfun3(TEST);
ERROR(2,3) = norm(GroundTruth - TEST);
disp('medium')

% Higher size.
tic
TEST = Heat3DMatEqnFD(f,100,20,.01);
TIME(3,3) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'11t');
TEST = chebfun3(TEST);
ERROR(3,3) = norm(GroundTruth - TEST);
disp('large')

% Higher size.
tic
TEST = Heat3DMatEqnFD(f,120,20,.01);
TIME(4,3) = toc;
TEST = equi2cheb(TEST(:,:,:,20),'11t');
TEST = chebfun3(TEST);
ERROR(4,3) = norm(GroundTruth - TEST);
disp('very large')