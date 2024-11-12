function result = evaluate(obj, fun, vertices)

[localNodes, Jacobi] = obj.localNodes(vertices);

temp = arrayfun(fun, localNodes(:,1), localNodes(:,2),'UniformOutput',false);
ele = cell2mat(temp(1));
pointRes = cell2mat(temp);

widenSize = size(ele,1);
weightMatrix = kron(obj.weight', eye(widenSize));


result = Jacobi * weightMatrix * pointRes;

end

