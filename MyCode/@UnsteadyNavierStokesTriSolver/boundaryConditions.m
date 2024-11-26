function boundaryConditions(obj,t)

for i = 1 : obj.boundaryInfo.numNodes
    
    
    % Dirichlet
    if obj.boundaryInfo.nodes(i,1) == -1
        index = obj.boundaryInfo.nodes(i,3);
        point = obj.uTrialElementInfo.P(index,:);
        obj.lhsMatrix(index,:) = 0;
        obj.lhsMatrix(index,index) = 1;
        obj.rhsVector(index,1) = obj.boundF1(point(1,1),point(1,2),t);
    end
    
    if obj.boundaryInfo.nodes(i,2) == -1
        index = obj.boundaryInfo.nodes(i,3);
        point = obj.uTrialElementInfo.P(index,:);
        rowSize = obj.uTestElementInfo.numPoints;
        colSize = obj.uTrialElementInfo.numPoints;
        obj.lhsMatrix(index + rowSize, :) = 0;
        obj.lhsMatrix(index + rowSize, index + colSize) = 1;
        obj.rhsVector(index + rowSize, 1) = obj.boundF2(point(1,1),point(1,2),t);
    end

end

end

