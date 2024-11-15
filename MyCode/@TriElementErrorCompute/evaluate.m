function diff = evaluate(obj,refFun,femSol,vertices,orderX,orderY)

[nodes,~] = obj.quadrature2D.localNodes(vertices);

fun = @(x,y) obj.trialFun.evaluate(x,y,vertices,orderX,orderY);

temp = arrayfun(fun, nodes(:,1), nodes(:,2),'UniformOutput',false);
pointRes = cell2mat(temp);
basis = reshape(pointRes,size(femSol,1),[]);

femGauss = basis' * femSol;

refGauss = arrayfun(refFun, nodes(:,1), nodes(:,2));

diff = refGauss - femGauss;

end

