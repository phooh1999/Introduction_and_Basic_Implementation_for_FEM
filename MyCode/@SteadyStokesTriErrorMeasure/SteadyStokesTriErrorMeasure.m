classdef SteadyStokesTriErrorMeasure < handle

    properties
        solution
        
        meshInfo
        uTrialElementInfo
        pTrialElementInfo
        
        uElementErrorCompute
        pElementErrorCompute
        
        f1
        f1_x
        f1_y
        
        f2
        f2_x
        f2_y
        
        p
        p_x
        p_y
    end
    
    methods
        function obj = SteadyStokesTriErrorMeasure(sol,mesh,uTrialEle,pTrialEle,uTrialFun,pTrialFun,quad,funRef1,funRef1_x,funRef1_y,funRef2,funRef2_x,funRef2_y,funRefp,funRefp_x,funRefp_y)
            obj.solution = sol;
            
            obj.meshInfo = mesh;
            obj.uTrialElementInfo = uTrialEle;
            obj.pTrialElementInfo = pTrialEle;
            
            obj.uElementErrorCompute = TriElementErrorCompute(uTrialFun,quad);
            obj.pElementErrorCompute = TriElementErrorCompute(pTrialFun,quad);
            
            obj.f1 = funRef1;
            obj.f1_x = funRef1_x;
            obj.f1_y = funRef1_y;
            
            obj.f2 = funRef2;
            obj.f2_x = funRef2_x;
            obj.f2_y = funRef2_y;
            
            obj.p = funRefp;
            obj.p_x = funRefp_x;
            obj.p_y = funRefp_y;
        end
        
        [errU,errP] = InfinityError(obj);
        [errU,errP] = L2Error(obj);
        [errU,errP] = H1Error(obj);
        
    end
end

