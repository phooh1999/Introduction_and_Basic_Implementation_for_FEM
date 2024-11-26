function [AN, bN] = generateNonlinearPart(obj,preSol)

P = obj.meshInfo.P;
T = obj.meshInfo.T;

% element size
uTrialEleSize = obj.uTrialElementInfo.numLocal;
uTestEleSize = obj.uTestElementInfo.numLocal;

% block matrix index
rowSize = obj.uTestElementInfo.numPoints;
colSize = obj.uTrialElementInfo.numPoints;
pRowSize = obj.pTestElementInfo.numPoints;
pColSize = obj.pTrialElementInfo.numPoints;

% index and value to generate sparse(lhs)
computeASize = (uTestEleSize*2)*(uTrialEleSize*2);
A_i = zeros(obj.numElements * computeASize, 1);
A_j = zeros(obj.numElements * computeASize, 1);
A_k = zeros(obj.numElements * computeASize, 1);

% index and value to generate sparse(rhs)
computebSize = uTestEleSize*2;
b_i = zeros(obj.numElements * computebSize,1);
b_j = ones(obj.numElements * computebSize,1);
b_k = zeros(obj.numElements * computebSize,1);

% iterate elements
for i = 1:obj.numElements
    
    vertices = P(T(i,:),:);
    solIndex = obj.uTrialElementInfo.T(i,:);
    femSolu1 = preSol(solIndex,1);
    femSolu2 = preSol(solIndex + colSize,1);
    
    u1 = @(x,y) obj.uuComputeMethod.getVal(x,y,0,0,femSolu1,vertices);
    u1_x = @(x,y) obj.uuComputeMethod.getVal(x,y,1,0,femSolu1,vertices);
    u1_y = @(x,y) obj.uuComputeMethod.getVal(x,y,0,1,femSolu1,vertices);
    
    u2 = @(x,y) obj.uuComputeMethod.getVal(x,y,0,0,femSolu2,vertices);
    u2_x = @(x,y) obj.uuComputeMethod.getVal(x,y,1,0,femSolu2,vertices);
    u2_y = @(x,y) obj.uuComputeMethod.getVal(x,y,0,1,femSolu2,vertices);
    
    % A index
    tmp1(:,1) = obj.uTestElementInfo.T(i,:);
    tmp2 = [tmp1;tmp1+rowSize];
    TA_test = kron(ones(2*uTestEleSize,1),tmp2);
    
    tmp3(:,1) = obj.uTrialElementInfo.T(i,:);
    tmp4 = [tmp3;tmp3+colSize];
    TA_trial = kron(tmp4,ones(2*uTrialEleSize,1));
    
    % element matrix
    AN1 = obj.uuComputeMethod.computeLhs(vertices,u1_x,0,0,0,0);
    AN2 = obj.uuComputeMethod.computeLhs(vertices,u1,1,0,0,0);
    AN3 = obj.uuComputeMethod.computeLhs(vertices,u2,0,1,0,0);
    AN4 = obj.uuComputeMethod.computeLhs(vertices,u1_y,0,0,0,0);
    AN5 = obj.uuComputeMethod.computeLhs(vertices,u2_x,0,0,0,0);
    AN6 = obj.uuComputeMethod.computeLhs(vertices,u2_y,0,0,0,0);

    
    AN = [AN1+AN2+AN3, AN4;
          AN5, AN6+AN2+AN3];
    
    AIndex = (i-1)*computeASize+1 : i*computeASize;
    
    % element matrix to total matrix
    A_i(AIndex) = TA_test;
    A_j(AIndex) = TA_trial;
    A_k(AIndex) = AN;
    
    % b index
    Tbu_test(:,1) = obj.uTestElementInfo.T(i,:);
    
    % element vector
    bN1 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) u1(x,y)*u1_x(x,y),0,0);
    bN2 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) u2(x,y)*u1_y(x,y),0,0);
    bN3 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) u1(x,y)*u2_x(x,y),0,0);
    bN4 = obj.uuComputeMethod.computeRhs(vertices,@(x,y) u2(x,y)*u2_y(x,y),0,0);
    b = [bN1+bN2;bN3+bN4];
    
    bIndex = (i-1)*computebSize+1 : i*computebSize;
    
    % element vector to total vector
    b_i(bIndex) = [Tbu_test;Tbu_test + rowSize];
    b_k(bIndex) = b;
    
end

AN = sparse(A_i,A_j,A_k,rowSize*2 + pRowSize,colSize*2 + pColSize);
bN = sparse(b_i,b_j,b_k,rowSize*2 + pRowSize,1);

end

