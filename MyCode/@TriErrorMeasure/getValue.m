function [diffX,diffY] = getValue(obj,funX,funY,index,orderX,orderY)

femVal1 = obj.solution.u1(obj.trialElementInfo.T(index,:),1);
femVal2 = obj.solution.u2(obj.trialElementInfo.T(index,:),1);

vertices = obj.meshInfo.P(obj.meshInfo.T(index,:),:);

[nodes,~] = obj.quadrature2D.localNodes(vertices);

fun = @(x,y) obj.trialFun.evaluate(x,y,vertices,orderX,orderY);

temp = arrayfun(fun, nodes(:,1), nodes(:,2),'UniformOutput',false);
pointRes = cell2mat(temp);
basis = reshape(pointRes,size(femVal1,1),[]);

% FEM value in each gauss points of ith element
femX = basis' * femVal1;
femY = basis' * femVal2;

% Exact value in each gauss points of ith element
refX = arrayfun(funX, nodes(:,1), nodes(:,2));
refY = arrayfun(funY, nodes(:,1), nodes(:,2));

diffX = refX - femX;
diffY = refY - femY;

end

