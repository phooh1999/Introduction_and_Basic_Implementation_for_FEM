classdef TriElementCompute < handle
    
    properties
        trialFun
        testFun
        
        quadrature2D
    end
    
    methods
        function obj = TriElementCompute(trial,test,gaussQuad)
            obj.trialFun = trial;
            obj.testFun = test;
            obj.quadrature2D = gaussQuad;
        end
        
        eleMatrix = computeLhs(obj,vertices,coe,trialOrderX,trialOrderY,testOrderX,testOrderY);
        eleVector = computeRhs(obj,vertices,coe,testOrderX,testOrderY);

        val = getVal(obj,x,y,orderX,orderY,femSol,vertices);
    end
end

