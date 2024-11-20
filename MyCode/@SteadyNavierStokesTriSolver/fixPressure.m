function fixPressure(obj,funRefp,index)

% fix p at one point in the domain Omega
point = obj.pTrialElementInfo.P(index,:);
rowSize = obj.uTestElementInfo.numPoints;
colSize = obj.uTrialElementInfo.numPoints;
obj.lhsMatrix(index + 2*rowSize, :) = 0;
obj.lhsMatrix(index + 2*rowSize, index + 2*colSize) = 1;
obj.rhsVector(index + 2*rowSize, 1) = funRefp(point(1,1),point(1,2));

end

