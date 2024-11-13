function err = L2Error(obj)

sum = 0;

for i = 1 : obj.meshInfo.numElements
    
    [diffX,diffY] = obj.getValue(obj.f1,obj.f2,i,0,0);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    [~,jac] = obj.quadrature2D.localNodes(vertices);
    
    sum = sum + jac * obj.quadrature2D.weight' * (diffX .* diffX);
    sum = sum + jac * obj.quadrature2D.weight' * (diffY .* diffY);
    
end

err = sqrt(sum);

end

