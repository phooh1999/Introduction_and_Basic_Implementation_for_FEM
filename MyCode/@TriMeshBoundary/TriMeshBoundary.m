classdef TriMeshBoundary < handle

    properties
        nodes
        edges
        numNodes
        numEdges
        
        partition
    end
    
    methods
        function obj = TriMeshBoundary(steps,type)
            generateBoundary(obj,steps,type);
        end

        setBoundaryType(obj,type,edge,coor);
        
    end
end

