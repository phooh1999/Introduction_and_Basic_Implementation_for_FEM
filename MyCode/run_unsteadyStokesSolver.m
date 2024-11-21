%% ---- Chapter 8: Finite elements for 2D unsteady Stokes equation ----
disp('----------------------------------------------')
formatSpec = '1/%2d, 1/%2d    %.4e    %.4e    %.4e\n';

%% Input
step = 1/32;
dt = 1/256;
theta = 1/2; % ODE solver type

% problem domain
problemDomain.x = [0,1];
problemDomain.y = [-0.25,0];

t0 = 0; t1 = 1;

% steps of the mesh partition
x_domain = problemDomain.x;
y_domain = problemDomain.y;
partitionSteps.x = (x_domain(2) - x_domain(1))/step;
partitionSteps.y = (y_domain(2) - y_domain(1))/step;

% basis function type
uTrialBasisType = 2;
uTestBasisType = 2;

pTrialBasisType = 1;
pTestBasisType = 1;

% number of gauss quadrature points
gaussPoints2D = 9;
gaussPointsLine = 4; % TODO: 边界条件处理中需要

% other functions required to solve the problem
funMu = @(x,y) 1;
funF1 = @(x,y,t) -2*pi*(x^2*y^2+exp(-y)) * sin(2*pi*t) + (-2*funMu(x,y)*x^2-2*funMu(x,y)*y^2-funMu(x,y)*exp(-y)+pi^2*cos(pi*x)*cos(2*pi*y))*cos(2*pi*t);
funF2 = @(x,y,t) -2*pi*(-2/3*x*y^3+2-pi*sin(pi*x))*sin(2*pi*t) + (4*funMu(x,y)*x*y - funMu(x,y)*pi^3*sin(pi*x) + 2*pi*(2-pi*sin(pi*x))*sin(2*pi*y))*cos(2*pi*t);

% boundary conditions
funBoundu1 = @bound1;
funBoundu2 = @bound2;

% initial value
iniu1 = @(x,y) x^2*y^2 + exp(-y);
iniu2 = @(x,y) -2/3*x*y^3 + 2 - pi*sin(pi*x);
inip = @(x,y) -(2-pi*sin(pi*x)) * cos(2*pi*y);

% reference solution
funRefu1_t = @(x,y,t) (x^2*y^2 + exp(-y)) * cos(2*pi*t);
funRefu1_x_t = @(x,y,t) 2*x*y^2 * cos(2*pi*t);
funRefu1_y_t = @(x,y,t) (2*y*x^2 - exp(-y)) * cos(2*pi*t);

funRefu2_t = @(x,y,t) (-2/3*x*y^3 + 2 - pi*sin(pi*x)) * cos(2*pi*t);
funRefu2_x_t = @(x,y,t) (-2/3*y^3 - pi^2*cos(pi*x)) * cos(2*pi*t);
funRefu2_y_t = @(x,y,t) -2*x*y^2 * cos(2*pi*t);

funRefp_t = @(x,y,t) -(2 - pi*sin(pi*x)) * cos(2*pi*y)  * cos(2*pi*t);
funRefp_x_t = @(x,y,t) pi^2*cos(pi*x)*cos(2*pi*y) * cos(2*pi*t);
funRefp_y_t = @(x,y,t) 2*pi*(2 - pi*sin(pi*x))*sin(2*pi*y) * cos(2*pi*t);

%% Solver
% Preprocessing
meshInfo = TriMeshGenerate(problemDomain, partitionSteps, 1);
uTrialElementInfo = TriMeshGenerate(problemDomain, partitionSteps, uTrialBasisType);
uTestElementInfo = TriMeshGenerate(problemDomain, partitionSteps, uTestBasisType);
pTrialElementInfo = TriMeshGenerate(problemDomain, partitionSteps, pTrialBasisType);
pTestElementInfo = TriMeshGenerate(problemDomain, partitionSteps, pTestBasisType);

boundary = TriMeshBoundary(partitionSteps, uTrialBasisType);

% Prepare Solver
uTrialFun = LocalTriBasisFunction(uTrialBasisType);
uTestFun = LocalTriBasisFunction(uTestBasisType);
pTrialFun = LocalTriBasisFunction(pTrialBasisType);
pTestFun = LocalTriBasisFunction(pTestBasisType);

quad = TriQuadrature(gaussPoints2D);
solver = UnsteadyStokesTriSolver(meshInfo,uTrialElementInfo,uTestElementInfo,pTrialElementInfo,pTestElementInfo,boundary,uTrialFun,uTestFun,pTrialFun,pTestFun,quad,funMu,funF1,funF2,funBoundu1,funBoundu2);

% Solve
solver.generateMatrix();

% solve ODE
x0 = solver.generateInitialVal(iniu1,iniu2,inip);
solver.solve(x0,t0,t1,dt,theta);

% Postprocessing
funRefu1 = @(x,y) funRefu1_t(x,y,t1);
funRefu1_x = @(x,y) funRefu1_x_t(x,y,t1);
funRefu1_y = @(x,y) funRefu1_y_t(x,y,t1);

funRefu2 = @(x,y) funRefu2_t(x,y,t1);
funRefu2_x = @(x,y) funRefu2_x_t(x,y,t1);
funRefu2_y = @(x,y) funRefu2_y_t(x,y,t1);

funRefp = @(x,y) funRefp_t(x,y,t1);
funRefp_x = @(x,y) funRefp_x_t(x,y,t1);
funRefp_y = @(x,y) funRefp_y_t(x,y,t1);

errMeasure = StokesTriErrorMeasure(solver.femSolution,meshInfo,uTrialElementInfo,pTrialElementInfo,uTrialFun,pTrialFun,quad,funRefu1,funRefu1_x,funRefu1_y,funRefu2,funRefu2_x,funRefu2_y,funRefp,funRefp_x,funRefp_y);

[uErrInfi,pErrInfi] = errMeasure.InfinityError();
[uErrL2,pErrL2] = errMeasure.L2Error();
[uErrH1,pErrH1] = errMeasure.H1Error();

%% Output
fprintf(formatSpec,1/step,1/dt,uErrInfi,uErrL2,uErrH1);
fprintf(formatSpec,1/step,1/dt,pErrInfi,pErrL2,pErrH1);

%% Appendix
function u1 = bound1(x,y,t)

if x == 0
    u1 = exp(-y) * cos(2*pi*t);
    return;
end

if x == 1
    u1 = (y^2 + exp(-y)) * cos(2*pi*t);
    return;
end

if y == -0.25
    u1 = (1/16*x^2 + exp(0.25)) * cos(2*pi*t);
    return;
end

if y == 0
    u1 = cos(2*pi*t);
    return;
end

end

function u2 = bound2(x,y,t)

if x == 0
    u2 = 2 * cos(2*pi*t);
    return;
end

if x == 1
    u2 = (-2/3*y^3 + 2) * cos(2*pi*t);
    return;
end

if y == -0.25
    u2 = (1/96*x + 2 - pi*sin(pi*x)) * cos(2*pi*t);
    return;
end

if y == 0
    u2 = (2 - pi*sin(pi*x)) * cos(2*pi*t);
    return;
end

end