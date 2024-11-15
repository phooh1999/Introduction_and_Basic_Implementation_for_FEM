function generateEquations(obj)

P = obj.meshInfo.P;
T = obj.meshInfo.T;

% element size
uTrialEleSize = obj.uTrialElementInfo.numLocal;
uTestEleSize = obj.uTestElementInfo.numLocal;
pTrialEleSize = obj.pTrialElementInfo.numLocal;
pTestEleSize = obj.pTestElementInfo.numLocal;

% block matrix index
rowSize = obj.uTestElementInfo.numPoints;
colSize = obj.uTrialElementInfo.numPoints;
pRowSize = obj.pTestElementInfo.numPoints;
pColSize = obj.pTrialElementInfo.numPoints;

% index and value to generate sparse(lhs)
computeASize = (uTestEleSize*2 + pTestEleSize)*(uTrialEleSize*2 + pTrialEleSize);
A_i = zeros(obj.numElements * computeASize, 1);
A_j = zeros(obj.numElements * computeASize, 1);
A_k = zeros(obj.numElements * computeASize, 1);

% index and value to generate sparse(rhs)
computebSize = uTestEleSize*2 + pTestEleSize;
b_i = zeros(obj.numElements * computebSize, 1);
b_j = ones(obj.numElements * computebSize,1);
b_k = zeros(obj.numElements * computebSize, 1);

% iterate elements
for i = 1:obj.numElements
    
    vertices = P(T(i,:),:);
    
    % A index
    tmp1(:,1) = obj.uTestElementInfo.T(i,:);
    tmp2(:,1) = obj.pTestElementInfo.T(i,:);
    tmp3 = [tmp1;tmp1+rowSize;tmp2+2*rowSize];
    TA_test = kron(ones(2*uTestEleSize + pTestEleSize,1),tmp3);
    
    tmp4(:,1) = obj.uTrialElementInfo.T(i,:);
    tmp5(:,1) = obj.pTrialElementInfo.T(i,:);
    tmp6 = [tmp4;tmp4+colSize;tmp5+2*colSize];
    TA_trial = kron(tmp6,ones(2*uTrialEleSize + pTrialEleSize,1));
    
    % element matrix
    A1 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,1,0,1,0);
    A2 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,0,1,0,1);
    A3 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,1,0,0,1);
    
    A5 = obj.puComputeMethod.computeLhs(vertices,@(x,y) -1,0,0,1,0);
    A6 = obj.puComputeMethod.computeLhs(vertices,@(x,y) -1,0,0,0,1);

    
    A = [2*A1+A2,      A3, A5;
             A3', A1+2*A2, A6;
             A5',     A6', zeros(pTestEleSize,pTrialEleSize)];
    
    AIndex = (i-1)*computeASize+1 : i*computeASize;
    
    % element matrix to total matrix
    A_i(AIndex) = TA_test;
    A_j(AIndex) = TA_trial;
    A_k(AIndex) = A;
    
    % b index
    Tbu_test(:,1) = obj.uTestElementInfo.T(i,:);
    Tbp_test(:,1) = obj.pTestElementInfo.T(i,:);
    b1 = obj.uuComputeMethod.computeRhs(vertices,obj.f1,0,0);
    b2 = obj.uuComputeMethod.computeRhs(vertices,obj.f2,0,0);
    b = [b1;b2;zeros(pTestEleSize,1)];
    
    bIndex = (i-1)*computebSize+1 : i*computebSize;
    
    % element vector to total vector
    b_i(bIndex) = [Tbu_test;Tbu_test + rowSize;Tbp_test + 2*rowSize];
    b_k(bIndex) = b;
    
end

obj.lhsMatrix = sparse(A_i,A_j,A_k,rowSize*2 + pRowSize,colSize*2 + pColSize);
obj.rhsVector = sparse(b_i,b_j,b_k,rowSize*2 + pRowSize,1);

end

