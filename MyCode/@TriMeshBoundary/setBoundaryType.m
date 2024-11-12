function setBoundaryType(obj,type,edge,coor)
% type: -1 -> Dirichlet boundary
%       -2 -> Stress boundary
%       -3 -> Robin boundary

% edge: 1 -> left
%       2 -> top
%       3 -> right
%       4 -> bottom

% coor: 1 -> x
%       2 -> y

nodeIndex = [];
edgeIndex = [];

xMeshStep = obj.partition(1);
yMeshStep = obj.partition(2);
xEleStep  = obj.partition(3);
yEleStep  = obj.partition(4);

if edge == 1
    nodeIndex = 1 : yEleStep;
    edgeIndex = 1 : yMeshStep;
elseif edge == 2
    nodeIndex = yEleStep + 1 : (xEleStep + yEleStep);
    edgeIndex = yMeshStep + 1 : xMeshStep + yMeshStep;
elseif edge == 3
    nodeIndex = (xEleStep + yEleStep + 1) : (xEleStep + 2*yEleStep);
    edgeIndex = xMeshStep + yMeshStep + 1 : xMeshStep + 2*yMeshStep;
elseif edge == 4
    nodeIndex = (xEleStep + 2*yEleStep + 1) : 2*(xEleStep + yEleStep);
    edgeIndex = xMeshStep + 2*yMeshStep + 1 : 2*(xMeshStep + yMeshStep);
else
    warning('Wrong TYPE of edge!');
end

obj.nodes(nodeIndex,coor) = type;
obj.edges(edgeIndex,coor) = type;

end

