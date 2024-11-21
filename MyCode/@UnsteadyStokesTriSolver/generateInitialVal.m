function x0 = generateInitialVal(obj,fu1,fu2,fp)

uSize = obj.uTrialElementInfo.numPoints;
pSize = obj.pTrialElementInfo.numPoints;
x0 = zeros(2*uSize + pSize,1);

uP = obj.uTrialElementInfo.P;
pP = obj.pTrialElementInfo.P;

for i = 1 : uSize
    
    x0(i,1) = fu1(uP(i,1),uP(i,2));
    x0(i+uSize,1) = fu2(uP(i,1),uP(i,2));
    
end

for j = 1 : pSize
    
    x0(j+2*uSize,1) = fp(pP(j,1),pP(j,2));
    
end 

x0 = sparse(x0);

end

