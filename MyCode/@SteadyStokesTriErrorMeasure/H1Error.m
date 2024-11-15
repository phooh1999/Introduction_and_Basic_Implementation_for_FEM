function [errU,errP] = H1Error(obj)

sumU = 0;
sumP = 0;

for i = 1 : obj.meshInfo.numElements
    
    femValu1 = obj.solution.u1(obj.uTrialElementInfo.T(i,:),1);
    femValu2 = obj.solution.u2(obj.uTrialElementInfo.T(i,:),1);
    femValp = obj.solution.p(obj.pTrialElementInfo.T(i,:),1);
    
    vertices = obj.meshInfo.P(obj.meshInfo.T(i,:),:);
    
    diffXX = obj.uElementErrorCompute.evaluate(obj.f1_x,femValu1,vertices,1,0);
    diffYX = obj.uElementErrorCompute.evaluate(obj.f2_x,femValu2,vertices,1,0);
    diffXY = obj.uElementErrorCompute.evaluate(obj.f1_y,femValu1,vertices,0,1);
    diffYY = obj.uElementErrorCompute.evaluate(obj.f2_y,femValu2,vertices,0,1);
    
    diffPX = obj.pElementErrorCompute.evaluate(obj.p_x,femValp,vertices,1,0);
    diffPY = obj.pElementErrorCompute.evaluate(obj.p_y,femValp,vertices,0,1);
    
    [~,jacU] = obj.uElementErrorCompute.quadrature2D.localNodes(vertices);
    [~,jacP] = obj.pElementErrorCompute.quadrature2D.localNodes(vertices);
    
    
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffXX .* diffXX);
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffYX .* diffYX);
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffXY .* diffXY);
    sumU = sumU + jacU * obj.uElementErrorCompute.quadrature2D.weight' * (diffYY .* diffYY);
    
    sumP = sumP + jacP * obj.pElementErrorCompute.quadrature2D.weight' * (diffPX .* diffPX);
    sumP = sumP + jacP * obj.pElementErrorCompute.quadrature2D.weight' * (diffPY .* diffPY);
    
end

errU = sqrt(sumU);
errP = sqrt(sumP);

end

