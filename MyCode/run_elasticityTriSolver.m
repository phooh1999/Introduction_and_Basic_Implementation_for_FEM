%% ---- Elasticity Triangle Solver ----
disp('----------------------------------------------')
disp('  h        Infi           L2            H1')
formatSpec = '1/%2d    %.4e    %.4e    %.4e\n';

tic
for step = [8 16 32 64]
%% Input
% problem domain
problemDomain.x = [0,1];
problemDomain.y = [0,1];

% steps of the mesh partition
partitionSteps.x = step;
partitionSteps.y = step;

% basis function type
trialBasisType = 2;
testBasisType = 2;

% number of gauss quadrature points
gaussPoints2D = 9;
gaussPointsLine = 4; % TODO: 边界条件处理中需要

% other functions required to solve the problem
funLambda = @(x,y) 1;
funMu = @(x,y) 2;
funF1 = @(x,y) -(funLambda(x,y)+2*funMu(x,y))*(-pi^2*sin(pi*x)*sin(pi*y))-(funLambda(x,y)+funMu(x,y))*((2*x-1)*(2*y-1))-funMu(x,y)*(-pi^2*sin(pi*x)*sin(pi*y));
funF2 = @(x,y) -(funLambda(x,y)+2*funMu(x,y))*(2*x*(x-1))-(funLambda(x,y)+funMu(x,y))*(pi^2*cos(pi*x)*cos(pi*y))-funMu(x,y)*(2*y*(y-1));

% boundary conditions
funBound1 = @(x,y) 0;
funBound2 = @(x,y) 0;

% reference solution
funRef1 = @(x,y) sin(pi*x)*sin(pi*y);
funRef1_x = @(x,y) pi*cos(pi*x)*sin(pi*y);
funRef1_y = @(x,y) pi*sin(pi*x)*cos(pi*y);

funRef2 = @(x,y) x*(x-1)*y*(y-1);
funRef2_x = @(x,y) (2*x-1)*y*(y-1);
funRef2_y = @(x,y) x*(x-1)*(2*y-1);

%% Solver
% Preprocessing
meshInfo = TriMeshGenerate(problemDomain, partitionSteps, 1);
trialElementInfo = TriMeshGenerate(problemDomain, partitionSteps, trialBasisType);
testElementInfo = TriMeshGenerate(problemDomain, partitionSteps, testBasisType);

boundary = TriMeshBoundary(partitionSteps, trialBasisType);

% Prepare Solver
trialFun = LocalTriBasisFunction(trialBasisType);
testFun = LocalTriBasisFunction(testBasisType);
quad = TriQuadrature(gaussPoints2D);
solver = ElasticityTriSolver(meshInfo,trialElementInfo,testElementInfo,boundary,trialFun,testFun,quad,funLambda,funMu,funF1,funF2,funBound1,funBound2);

% Solve
solver.generateEquations();
solver.boundaryConditions();
solver.solve();

% Postprocessing
errMeasure = ElasticityTriErrorMeasure(solver.femSolution,meshInfo,trialElementInfo,trialFun,quad,funRef1,funRef1_x,funRef1_y,funRef2,funRef2_x,funRef2_y);

errInfi = errMeasure.InfinityError();
errL2 = errMeasure.L2Error();
errH1 = errMeasure.H1Error();

%% Output
fprintf(formatSpec,partitionSteps.x,errInfi,errL2,errH1)
end
toc
disp('----------------------------------------------')
