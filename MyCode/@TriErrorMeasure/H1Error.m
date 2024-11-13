function err = H1Error(obj)

sum = 0;

for i = 1 : obj.meshInfo.numElements
    
    [diffXX,diffYX] = obj.getValue(obj.f1_x,obj.f2_x,i,1,0);
    [diffXY,diffYY] = obj.getValue(obj.f1_y,obj.f2_y,i,0,1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    [~,jac] = obj.quadrature2D.localNodes(vertices);
    
    sum = sum + jac * obj.quadrature2D.weight' * (diffXX .* diffXX);
    sum = sum + jac * obj.quadrature2D.weight' * (diffYX .* diffYX);
    sum = sum + jac * obj.quadrature2D.weight' * (diffXY .* diffXY);
    sum = sum + jac * obj.quadrature2D.weight' * (diffYY .* diffYY);
    
end

err = sqrt(sum);

end

