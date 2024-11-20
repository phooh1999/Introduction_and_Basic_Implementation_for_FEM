function [errU,errP] = InfinityError(obj)

maxUX = 0;
maxUY = 0;
maxP = 0;

for i = 1 : obj.meshInfo.numElements
    
    femValu1 = obj.solution.u1(obj.uTrialElementInfo.T(i,:),1);
    femValu2 = obj.solution.u2(obj.uTrialElementInfo.T(i,:),1);
    femValp = obj.solution.p(obj.pTrialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffUX = obj.uElementErrorCompute.evaluate(obj.f1,femValu1,vertices,0,0);
    diffUY = obj.uElementErrorCompute.evaluate(obj.f2,femValu2,vertices,0,0);
    diffP = obj.pElementErrorCompute.evaluate(obj.p,femValp,vertices,0,0);
    
    tempUX = max(abs(diffUX));
    tempUY = max(abs(diffUY));
    tempP = max(abs(diffP));
    
    if tempUX > maxUX
        maxUX = tempUX;
    end
    
    if tempUY > maxUY
        maxUY = tempUY;
    end
    
    if tempP > maxP
        maxP = tempP;
    end
    
end

errU = max([maxUX,maxUY]);
errP = maxP;

end

