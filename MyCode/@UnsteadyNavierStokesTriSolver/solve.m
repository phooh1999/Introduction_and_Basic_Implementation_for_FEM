function solve(obj,x0,t0,t1,dt,theta,maxStep,rtol,atol)

N = (t1-t0)/dt;

lhsLinearMat = obj.MMatrix ./dt + theta .* obj.AMatrix;
rhsLinearMat = obj.MMatrix ./dt - (1-theta) .* obj.AMatrix;

x_m = x0;
b_m = obj.generateVector(t0);
[~,bN_m] = obj.generateNonlinearPart(x_m);

% fix pressure at one point
funRefp = @(x,y,t) -(2 - pi*sin(pi*x)) * cos(2*pi*y)  * cos(2*pi*t);

% solve ODEs
for i = 1:N
    
    time = t0 + i*dt;
    
    b_m1 = obj.generateVector(time);
    rhsLinearVec = theta*b_m1 + (1-theta)*b_m - (1-theta)*bN_m + rhsLinearMat*x_m;
    
    x_m1 = x_m;
    
    % Newton's iteration
    for j = 1:maxStep
    
        
        [AN_m1,bN_m1] = obj.generateNonlinearPart(x_m1);

        obj.lhsMatrix = lhsLinearMat + theta*AN_m1;
        obj.rhsVector = rhsLinearVec + theta*bN_m1;

        obj.boundaryConditions(time);

        % for Dirichlet boundary conditions
        % fix pressure at one point in the domain Omega
        obj.fixPressure(funRefp,1,time);

        sol = obj.lhsMatrix\obj.rhsVector;

        % error estimation
        Wq = rtol * abs(sol) + atol;
        norm = max(abs(sol - x_m1) ./ Wq);

        x_m1 = sol;

        if norm < 1e-4
            break
        end
    
    end
    
    x_m = x_m1;
    b_m = b_m1;
    bN_m = bN_m1;
    
end

uColSize = obj.uTrialElementInfo.numPoints;
pColSize = obj.pTrialElementInfo.numPoints;
solution.u1 = x_m(1:uColSize);
solution.u2 = x_m(uColSize+1:2*uColSize);
solution.p = x_m(2*uColSize+1:2*uColSize+pColSize);
obj.femSolution = solution;

end
