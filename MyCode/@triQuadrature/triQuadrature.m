classdef TriQuadrature < handle
    
    properties
        lroots
        weight
    end
    
    methods
        function obj = TriQuadrature(points)
            tablesIntegrateTriangle(obj, points);
        end
        
        result = evaluate(obj, fun, vertices)
        
    end
end

