%% ---- Chapter 7: Finite elements for 2D steady Navier-Stokes equation ----
disp('----------------------------------------------')
formatSpec = '1/%2d    %.4e    %.4e    %.4e\n';

uErrInfi = zeros(4,1); uErrL2 = zeros(4,1); uErrH1 = zeros(4,1);
pErrInfi = zeros(4,1); pErrL2 = zeros(4,1); pErrH1 = zeros(4,1);

for i = 1:4
%% Input
step = 1/(2^(i+2));
% problem domain
problemDomain.x = [0,1];
problemDomain.y = [-0.25,0];

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
funF1 = @(x,y) -2*funMu(x,y)*x^2 - 2*funMu(x,y)*y^2 - funMu(x,y)*exp(-y) + pi^2*cos(pi*x)*cos(2*pi*y) + 2*x*y^2*(x^2*y^2+exp(-y)) + (-2*x*y^3/3+2-pi*sin(pi*x))*(2*x^2*y-exp(-y));
funF2 = @(x,y) 4*funMu(x,y)*x*y - funMu(x,y)*pi^3*sin(pi*x) + 2*pi*(2-pi*sin(pi*x))*sin(2*pi*y) + (x^2*y^2+exp(-y))*(-2*y^3/3-pi^2*cos(pi*x)) + (-2*x*y^3/3+2-pi*sin(pi*x))*(-2*x*y^2);

% boundary conditions
funBoundu1 = @bound1;
funBoundu2 = @bound2;

% reference solution
funRefu1 = @(x,y) x^2*y^2 + exp(-y);
funRefu1_x = @(x,y) 2*x*y^2;
funRefu1_y = @(x,y) 2*y*x^2 - exp(-y);

funRefu2 = @(x,y) -2/3*x*y^3 + 2 - pi*sin(pi*x);
funRefu2_x = @(x,y) -2/3*y^3 - pi^2*cos(pi*x);
funRefu2_y = @(x,y) -2*x*y^2;

funRefp = @(x,y) -(2 - pi*sin(pi*x)) * cos(2*pi*y);
funRefp_x = @(x,y) pi^2*cos(pi*x)*cos(2*pi*y);
funRefp_y = @(x,y) 2*pi*(2 - pi*sin(pi*x))*sin(2*pi*y);

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
solver = SteadyNavierStokesTriSolver(meshInfo,uTrialElementInfo,uTestElementInfo,pTrialElementInfo,pTestElementInfo,boundary,uTrialFun,uTestFun,pTrialFun,pTestFun,quad,funMu,funF1,funF2,funBoundu1,funBoundu2);

% initial solution guess
uColSize = solver.uTrialElementInfo.numPoints;
pColSize = solver.pTrialElementInfo.numPoints;
iniSol = sparse(2*uColSize+pColSize,1);

% Solve
solver.generateLinearPart();
solver.solve(iniSol,25,1e-3,1e-6);

% Postprocessing
errMeasure = StokesTriErrorMeasure(solver.femSolution,meshInfo,uTrialElementInfo,pTrialElementInfo,uTrialFun,pTrialFun,quad,funRefu1,funRefu1_x,funRefu1_y,funRefu2,funRefu2_x,funRefu2_y,funRefp,funRefp_x,funRefp_y);

[uErrInfi(i),pErrInfi(i)] = errMeasure.InfinityError();
[uErrL2(i),pErrL2(i)] = errMeasure.L2Error();
[uErrH1(i),pErrH1(i)] = errMeasure.H1Error();
end

%% Output
for i = 1:4
fprintf(formatSpec,2^(i+2),uErrInfi(i),uErrL2(i),uErrH1(i));
end
disp('----------------------------------------------')
for i = 1:4
fprintf(formatSpec,2^(i+2),pErrInfi(i),pErrL2(i),pErrH1(i));
end
disp('----------------------------------------------')

%% Appendix
function u1 = bound1(x,y)

if x == 0
    u1 = exp(-y);
    return;
end

if x == 1
    u1 = y^2 + exp(-y);
    return;
end

if y == -0.25
    u1 = 1/16*x^2 + exp(0.25);
    return;
end

if y == 0
    u1 = 1;
    return;
end

end

function u2 = bound2(x,y)

if x == 0
    u2 = 2;
    return;
end

if x == 1
    u2 = -2/3*y^3 + 2;
    return;
end

if y == -0.25
    u2 = 1/96*x + 2 - pi*sin(pi*x);
    return;
end

if y == 0
    u2 = 2 - pi*sin(pi*x);
    return;
end

end