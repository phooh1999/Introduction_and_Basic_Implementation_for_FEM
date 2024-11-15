function err = L2Error(obj)

sum = 0;

for i = 1 : obj.meshInfo.numElements
    
    femVal1 = obj.solution.u1(obj.trialElementInfo.T(i,:),1);
    femVal2 = obj.solution.u2(obj.trialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffX = obj.elementErrorCompute.evaluate(obj.f1,femVal1,vertices,0,0);
    diffY = obj.elementErrorCompute.evaluate(obj.f2,femVal2,vertices,0,0);
    
    [~,jac] = obj.elementErrorCompute.quadrature2D.localNodes(vertices);
    
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffX .* diffX);
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffY .* diffY);
    
end

err = sqrt(sum);

end

