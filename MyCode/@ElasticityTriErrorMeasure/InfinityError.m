function err = InfinityError(obj)

maxX = 0;
maxY = 0;

for i = 1 : obj.meshInfo.numElements
    
    femVal1 = obj.solution.u1(obj.trialElementInfo.T(i,:),1);
    femVal2 = obj.solution.u2(obj.trialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffX = obj.elementErrorCompute.evaluate(obj.f1,femVal1,vertices,0,0);
    diffY = obj.elementErrorCompute.evaluate(obj.f2,femVal2,vertices,0,0);
    
    tempX = max(abs(diffX));
    tempY = max(abs(diffY));
    
    if tempX > maxX
        maxX = tempX;
    end
    
    if tempY > maxY
        maxY = tempY;
    end
    
end

err = max([maxX,maxY]);

end

