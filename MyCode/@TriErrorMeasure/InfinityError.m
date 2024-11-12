function err = InfinityError(obj)

maxX = 0;
maxY = 0;

for i = 1 : obj.meshInfo.numElements
    
    [diffX,diffY] = obj.getValue(obj.f1,obj.f2,i,0,0);
    
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

