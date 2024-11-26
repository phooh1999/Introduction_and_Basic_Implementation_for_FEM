function bVec = generateVector(obj,t)


P = obj.meshInfo.P;
T = obj.meshInfo.T;

% element size
uTestEleSize = obj.uTestElementInfo.numLocal;

% block matrix index
rowSize = obj.uTestElementInfo.numPoints;
pRowSize = obj.pTestElementInfo.numPoints;

% index and value to generate sparse(rhs)
computebSize = uTestEleSize*2;
b_i = zeros(obj.numElements * computebSize, 1);
b_j = ones(obj.numElements * computebSize,1);
b_k = zeros(obj.numElements * computebSize, 1);

% iterate elements
for i = 1:obj.numElements
    
    vertices = P(T(i,:),:);
    
    % b index
    Tbu_test(:,1) = obj.uTestElementInfo.T(i,:);
    b1 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) obj.f1(x,y,t),0,0);
    b2 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) obj.f2(x,y,t),0,0);
    b = [b1;b2];
    
    bIndex = (i-1)*computebSize+1 : i*computebSize;
    
    % element vector to total vector
    b_i(bIndex) = [Tbu_test;Tbu_test + rowSize];
    b_k(bIndex) = b;
    
end

bVec = sparse(b_i,b_j,b_k,rowSize*2 + pRowSize,1);

end

