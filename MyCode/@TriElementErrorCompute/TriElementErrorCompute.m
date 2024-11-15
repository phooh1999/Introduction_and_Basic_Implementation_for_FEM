classdef TriElementErrorCompute < handle
    
    properties
        trialFun

        quadrature2D
    end
    
    methods
        function obj = TriElementErrorCompute(trial,gaussQuad)
            obj.trialFun = trial;
            obj.quadrature2D = gaussQuad;
        end
        
        diff = evaluate(obj,refFun,femSol,vertices,orderX,orderY);

    end
end

