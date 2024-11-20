function solve(obj,iniSol,maxStep,rtol,atol)

% initial guess
preSol = iniSol;

% fix pressure
funRefp = @(x,y) -(2 - pi*sin(pi*x)) * cos(2*pi*y);

% Newton's iteration
for i = 1:maxStep
    
    [AN,bN] = generateNonlinearPart(obj,preSol);
    
    obj.lhsMatrix = obj.lhsLinearMatrix + AN;
    obj.rhsVector = obj.rhsLinearVector + bN;
    
    obj.boundaryConditions();
    
    % for Dirichlet boundary conditions
    % fix pressure at one point in the domain Omega
    obj.fixPressure(funRefp,1);
    
    sol = obj.lhsMatrix\obj.rhsVector;
    
    % error estimation
    Wq = rtol * abs(sol) + atol;
    norm = max(abs(sol - preSol) ./ Wq);
    
    preSol = sol;
    
    if norm < 1e-4
        break
    end
    
end

uColSize = obj.uTrialElementInfo.numPoints;
pColSize = obj.pTrialElementInfo.numPoints;
solution.u1 = preSol(1:uColSize);
solution.u2 = preSol(uColSize+1:2*uColSize);
solution.p = preSol(2*uColSize+1:2*uColSize+pColSize);
obj.femSolution = solution;

end

