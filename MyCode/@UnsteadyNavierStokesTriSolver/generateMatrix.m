function generateMatrix(obj)

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

% index and value to generate the sparse A matrix
computeASize = (uTestEleSize*2 + pTestEleSize)*(uTrialEleSize*2 + pTrialEleSize);
A_i = zeros(obj.numElements * computeASize, 1);
A_j = zeros(obj.numElements * computeASize, 1);
A_k = zeros(obj.numElements * computeASize, 1);

% index and value to generate the sparse mass matrix
computeMSize = (uTestEleSize)*(uTrialEleSize)*2;
M_i = zeros(obj.numElements * computeMSize, 1);
M_j = zeros(obj.numElements * computeMSize, 1);
M_k = zeros(obj.numElements * computeMSize, 1);

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
    
    % element A matrix
    A1 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,1,0,1,0);
    A2 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,0,1,0,1);
    A3 = obj.uuComputeMethod.computeLhs(vertices,obj.fMu,1,0,0,1);
    
    A5 = obj.puComputeMethod.computeLhs(vertices,@(x,y) -1,0,0,1,0);
    A6 = obj.puComputeMethod.computeLhs(vertices,@(x,y) -1,0,0,0,1);

    
    A = [2*A1+A2,      A3, A5;
             A3', A1+2*A2, A6;
             A5',     A6', zeros(pTestEleSize,pTrialEleSize)];
    
    AIndex = (i-1)*computeASize+1 : i*computeASize;
    
    % element matrix to total A matrix
    A_i(AIndex) = TA_test;
    A_j(AIndex) = TA_trial;
    A_k(AIndex) = A;
    
    % M index
    tmpM1(:,1) = obj.uTestElementInfo.T(i,:);
    TM_test = kron(ones(uTestEleSize,1),tmpM1);
    
    tmpM2(:,1) = obj.uTrialElementInfo.T(i,:);
    TM_trial = kron(tmpM2,ones(uTrialEleSize,1));
    
    % element M matrix
    Me = obj.uuComputeMethod.computeLhs(vertices,@(x,y) 1,0,0,0,0);
    MIndex = (i-1)*computeMSize+1 : i*computeMSize;
    
    % element M matrix to total M matrix
    M_i(MIndex) = [TM_test;TM_test+rowSize];
    M_j(MIndex) = [TM_trial;TM_trial+colSize];
    M_k(MIndex) = [Me,Me];
    
end

obj.AMatrix = sparse(A_i,A_j,A_k,rowSize*2 + pRowSize,colSize*2 + pColSize);
obj.MMatrix = sparse(M_i,M_j,M_k,rowSize*2 + pRowSize,colSize*2 + pColSize);

end

