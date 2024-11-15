function err = H1Error(obj)

sum = 0;

for i = 1 : obj.meshInfo.numElements
    
    femVal1 = obj.solution.u1(obj.trialElementInfo.T(i,:),1);
    femVal2 = obj.solution.u2(obj.trialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffXX = obj.elementErrorCompute.evaluate(obj.f1_x,femVal1,vertices,1,0);
    diffYX = obj.elementErrorCompute.evaluate(obj.f2_x,femVal2,vertices,1,0);
    diffXY = obj.elementErrorCompute.evaluate(obj.f1_y,femVal1,vertices,0,1);
    diffYY = obj.elementErrorCompute.evaluate(obj.f2_y,femVal2,vertices,0,1);
    
    [~,jac] = obj.elementErrorCompute.quadrature2D.localNodes(vertices);
    
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffXX .* diffXX);
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffYX .* diffYX);
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffXY .* diffXY);
    sum = sum + jac * obj.elementErrorCompute.quadrature2D.weight' * (diffYY .* diffYY);
    
end

err = sqrt(sum);

end

