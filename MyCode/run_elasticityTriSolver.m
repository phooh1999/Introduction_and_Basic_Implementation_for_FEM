%% ---- Elasticity Triangle Solver ----

%% Input
% problem domain
problemDomain.x = [0,1];
problemDomain.y = [0,1];

% steps of the mesh partition
partitionSteps.x = 8;
partitionSteps.y = 8;

% basis function type
trialBasisType = 2;
testBasisType = 2;

% number of gauss quadrature points
gaussPoints2D = 9;
gaussPointsLine = 4; % TODO: 这个是边界条件处理中需要的，需要留一个处理接口

% other functions required to solve the problem
funLambda = @(x,y) 1;
funMu = @(x,y) 2;
funF1 = @(x,y) -(funLambda(x,y)+2*funMu(x,y))*(-pi^2*sin(pi*x)*sin(pi*y))-(funLambda(x,y)+funMu(x,y))*((2*x-1)*(2*y-1))-funMu(x,y)*(-pi^2*sin(pi*x)*sin(pi*y));
funF2 = @(x,y) -(funLambda(x,y)+2*funMu(x,y))*(2*x*(x-1))-(funLambda(x,y)+funMu(x,y))*(pi^2*cos(pi*x)*cos(pi*y))-funMu(x,y)*(2*y*(y-1));
% boundary conditions


% reference solution


%% Solver
% Preprocessing
meshInfo = TriMeshGenerate(problemDomain, partitionSteps, 1);
trialElementInfo = TriMeshGenerate(problemDomain, partitionSteps, trialBasisType);
testElementInfo = TriMeshGenerate(problemDomain, partitionSteps, testBasisType);

% Prepare Solver
trialFun = LocalTriBasisFunction(trialBasisType);
testFun = LocalTriBasisFunction(testBasisType);
quad = TriQuadrature(gaussPoints2D);
solver = ElasticityTriSolver(meshInfo,trialElementInfo,testElementInfo,trialFun,testFun,quad,funLambda,funMu,funF1,funF2);

% Solve
solver.generateEquation();


% Postprocessing




%% Output


