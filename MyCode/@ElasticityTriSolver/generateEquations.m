function generateEquations(obj)

P = obj.meshInfo.P;
T = obj.meshInfo.T;

% element size
trialEleSize = obj.trialElementInfo.numLocal;
testEleSize = obj.testElementInfo.numLocal;

% block matrix index
rowSize = obj.testElementInfo.numPoints;
colSize = obj.trialElementInfo.numPoints;

% index and value to generate sparse(lhs)
computeASize = testEleSize * 2 * trialEleSize * 2;
A_i = zeros(obj.numElements * computeASize, 1);
A_j = zeros(obj.numElements * computeASize, 1);
A_k = zeros(obj.numElements * computeASize, 1);

% index and value to generate sparse(rhs)
computebSize = testEleSize * 2;
b_i = zeros(obj.numElements * computebSize, 1);
b_j = ones(obj.numElements * computebSize,1);
b_k = zeros(obj.numElements * computebSize, 1);

% iterate elements
for i = 1:obj.numElements
    
    vertices = P(T(i,:),:);
    
    % A index
    tmp1(:,1) = obj.testElementInfo.T(i,:);
    tmp2 = [tmp1;tmp1+rowSize];
    TA_test = kron(ones(2*testEleSize,1),tmp2);
    
    tmp3(:,1) = obj.trialElementInfo.T(i,:);
    tmp4 = [tmp3;tmp3+colSize];
    TA_trial = kron(tmp4,ones(2*trialEleSize,1));
    
    % element matrix
    A1 = obj.elementComputeMethod.computeLhs(vertices,obj.fLambda,1,0,1,0);
    A2 = obj.elementComputeMethod.computeLhs(vertices,obj.fMu,1,0,1,0);
    A3 = obj.elementComputeMethod.computeLhs(vertices,obj.fMu,0,1,0,1);
    A4 = obj.elementComputeMethod.computeLhs(vertices,obj.fLambda,0,1,1,0);
    A5 = obj.elementComputeMethod.computeLhs(vertices,obj.fMu,1,0,0,1);
    A6 = obj.elementComputeMethod.computeLhs(vertices,obj.fLambda,1,0,0,1);
    A7 = obj.elementComputeMethod.computeLhs(vertices,obj.fMu,0,1,1,0);
    A8 = obj.elementComputeMethod.computeLhs(vertices,obj.fLambda,0,1,0,1);
    
    A = [A1+2*A2+A3,     A4+A5;
              A6+A7,A8+2*A3+A2];
    
    AIndex = (i-1)*computeASize+1 : i*computeASize;
    
    % element matrix to total matrix
    A_i(AIndex) = TA_test;
    A_j(AIndex) = TA_trial;
    A_k(AIndex) = A;
    
    % b index
    Tb_test(:,1) = obj.testElementInfo.T(i,:);
    b1 = obj.elementComputeMethod.computeRhs(vertices,obj.f1,0,0);
    b2 = obj.elementComputeMethod.computeRhs(vertices,obj.f2,0,0);
    b = [b1;b2];
    
    bIndex = (i-1)*computebSize+1 : i*computebSize;
    
    % element vector to total vector
    b_i(bIndex) = [Tb_test;Tb_test + rowSize];
    b_k(bIndex) = b;
    
end

obj.lhsMatrix = sparse(A_i,A_j,A_k,rowSize*2,colSize*2);
obj.rhsVector = sparse(b_i,b_j,b_k,rowSize*2,1);

end

