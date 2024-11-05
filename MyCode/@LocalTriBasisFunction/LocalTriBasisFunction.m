classdef LocalTriBasisFunction < handle
    
    properties
        numBasis
        refBasisFun
    end
    
    methods
        function obj = LocalTriBasisFunction(type)
            refBasisLib(obj,type);
        end
        
        res = evaluate(obj,x,y,vertices,orderX,orderY);

    end
end

