function result = evaluate(obj, fun, vertices)
%EVALUATE 高斯积分求解函数
%
x1=vertices(1,1);   y1=vertices(1,2);
x2=vertices(2,1);   y2=vertices(2,2);
x3=vertices(3,1);   y3=vertices(3,2);
            
Jacobi=0.5 * abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

localNodes(:,1) = x1+(x2-x1)*obj.lroots(:,1)+(x3-x1)*obj.lroots(:,2);
localNodes(:,2) = y1+(y2-y1)*obj.lroots(:,1)+(y3-y1)*obj.lroots(:,2);

temp = arrayfun(fun, localNodes(:,1), localNodes(:,2),'UniformOutput',false);
ele = cell2mat(temp(1));
pointRes = cell2mat(temp);

widenSize = size(ele,1);
weightMatrix = kron(obj.weight', eye(widenSize));


result = Jacobi * weightMatrix * pointRes;

end

