classdef ElasticityTriErrorMeasure < handle

    properties
        solution
        
        meshInfo
        trialElementInfo
        
        elementErrorCompute
        
        f1
        f1_x
        f1_y
        
        f2
        f2_x
        f2_y
    end
    
    methods
        function obj = ElasticityTriErrorMeasure(sol,mesh,trialEle,trialFun,quad,funRef1,funRef1_x,funRef1_y,funRef2,funRef2_x,funRef2_y)
            obj.solution = sol;
            obj.meshInfo = mesh;
            obj.trialElementInfo = trialEle;
            
            obj.elementErrorCompute = TriElementErrorCompute(trialFun,quad);
            
            obj.f1 = funRef1;
            obj.f1_x = funRef1_x;
            obj.f1_y = funRef1_y;
            obj.f2 = funRef2;
            obj.f2_x = funRef2_x;
            obj.f2_y = funRef2_y;
        end
        
        err = InfinityError(obj);
        err = L2Error(obj);
        err = H1Error(obj);
        
    end
end

