classdef TriMeshGenerate < handle

    properties
        P
        T
        numPoints
        numElements
        numLocal
    end
    
    methods
        function obj = TriMeshGenerate(domain, steps, type)
            generatePT(obj, domain, steps, type);
            obj.numPoints = size(obj.P,1);
            obj.numElements = size(obj.T,1);
            obj.numLocal = size(obj.T,2);
        end
    end
end

