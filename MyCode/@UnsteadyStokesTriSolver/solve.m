function solve(obj,x0,t0,t1,dt,theta)

N = (t1-t0)/dt;

obj.lhsMatrix = obj.MMatrix ./dt + theta .* obj.AMatrix;
A_tob = obj.MMatrix ./dt - (1-theta) .* obj.AMatrix;

x_tm = x0;
b_tm = obj.generateVector(t0);

% fix pressure at one point
funRefp = @(x,y,t) -(2 - pi*sin(pi*x)) * cos(2*pi*y)  * cos(2*pi*t);

for i = 1:N
    
    time = t0 + i*dt;
    
    b_tm1 = obj.generateVector(time);
    
    obj.rhsVector = theta * b_tm1 + (1-theta) * b_tm + A_tob * x_tm;
    
    %boundary conditions
    obj.boundaryConditions(time);
    obj.fixPressure(funRefp,1,time);
    
    x_tm = obj.lhsMatrix\obj.rhsVector;
    b_tm = b_tm1;
    
end

uColSize = obj.uTrialElementInfo.numPoints;
pColSize = obj.pTrialElementInfo.numPoints;
solution.u1 = x_tm(1:uColSize);
solution.u2 = x_tm(uColSize+1:2*uColSize);
solution.p = x_tm(2*uColSize+1:2*uColSize+pColSize);
obj.femSolution = solution;

end
