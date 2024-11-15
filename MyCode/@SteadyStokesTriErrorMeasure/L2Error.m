function [errU,errP] = L2Error(obj)

sumU = 0;
sumP = 0;

for i = 1 : obj.meshInfo.numElements
    
    femValu1 = obj.solution.u1(obj.uTrialElementInfo.T(i,:),1);
    femValu2 = obj.solution.u2(obj.uTrialElementInfo.T(i,:),1);
    femValp = obj.solution.p(obj.pTrialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffUX = obj.uElementErrorCompute.evaluate(obj.f1,femValu1,vertices,0,0);
    diffUY = obj.uElementErrorCompute.evaluate(obj.f2,femValu2,vertices,0,0);
    diffP = obj.pElementErrorCompute.evaluate(obj.p,femValp,vertices,0,0);
    
    [~,jacU] = obj.uElementErrorCompute.quadrature2D.localNodes(vertices);
    [~,jacP] = obj.pElementErrorCompute.quadrature2D.localNodes(vertices);
    
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffUX .* diffUX);
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffUY .* diffUY);
    sumP = sumP + jacP * obj.pElementErrorCompute.quadrature2D.weight' * (diffP .* diffP);
end

errU = sqrt(sumU);
errP = sqrt(sumP);

end

