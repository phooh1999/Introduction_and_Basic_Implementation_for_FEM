function boundaryConditions(obj)

for i = 1 : obj.boundaryInfo.numNodes
    
    
    % Dirichlet
    if obj.boundaryInfo.nodes(i,1) == -1
        index = obj.boundaryInfo.nodes(i,3);
        obj.lhsMatrix(index,:) = 0;
        obj.lhsMatrix(index,index) = 1;
        obj.rhsVector(index,1) = obj.boundF1(obj.testElementInfo.P(index,:));
    end
    
    if obj.boundaryInfo.nodes(i,2) == -1
        index = obj.boundaryInfo.nodes(i,3);
        colSize = obj.trialElementInfo.numPoints;
        obj.lhsMatrix(index + colSize, :) = 0;
        obj.lhsMatrix(index + colSize, index + colSize) = 1;
        obj.rhsVector(index + colSize, 1) = obj.boundF2(obj.testElementInfo.P(index,:));
    end

end

end

