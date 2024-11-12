classdef TriErrorMeasure < handle

    properties
        solution
        
        meshInfo
        trialElementInfo
        
        trialFun
        quadrature2D
        
        f1
        f1_x
        f1_y
        
        f2
        f2_x
        f2_y
    end
    
    methods
        function obj = TriErrorMeasure(sol,mesh,trialEle,trialFun,quad,funRef1,funRef1_x,funRef1_y,funRef2,funRef2_x,funRef2_y)
            obj.solution = sol;
            obj.meshInfo = mesh;
            obj.trialElementInfo = trialEle;
            obj.trialFun = trialFun;
            obj.quadrature2D = quad;
            obj.f1 = funRef1;
            obj.f1_x = funRef1_x;
            obj.f1_y = funRef1_y;
            obj.f2 = funRef2;
            obj.f2_x = funRef2_x;
            obj.f2_y = funRef2_y;
        end
        
        [valX,valY] = getValue(obj,funX,funY,index,orderX,orderY);
        err = InfinityError(obj);
        err = L2Error(obj);
        err = H1Error(obj);
        
    end
end

