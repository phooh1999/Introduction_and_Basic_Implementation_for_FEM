function generateBoundary(obj,steps,type)
% Boundary Nodes
% boundary nodes(1,k): specifiy the type of the kth boundary node for the normal direction(or u1).
% boundary nodes(1,k)=-1: Dirichlet boundary node in the normal direction(or u1);
% boundary nodes(1,k)=-2: Stress boundary node in the normal direction(or u1);
% boundary nodes(1,k)=-3: Robin boundary node in the normal direction(or u1). 
% boundary nodes(2,k): specifiy the type of the kth boundary node for the tangential direction(or u2).
% boundary nodes(2,k)=-1: Dirichlet boundary node in the tangential direction(or u2);
% boundary nodes(2,k)=-2: Stress boundary node in the tangential direction(or u2);
% boundary nodes(2,k)=-3: Robin boundary node in the tangential direction(or u2).
% The intersection node between Dirichlet boundary and other boundaries is a Dirichlet boundary node.
% boundary_nodes(3,k): global index of the kth boundary node among all nodes of FE. 
%                     That is, the index of FE is used here.
% boundary_nodes(4:5,k): the unit outer normal vector at the kth boundary node.
% boundary_nodes(6:7,k): the unit tangential vector at the kth boundary node. 
%                       It's counterclockwise to go from the normal vector to the tangential vector.

% Boundary Edges
% boundary_edges(1,k): specifiy the type of the kth boundary edge in normal direction.
% boundary_edges(1,k)=-1: Dirichlet boundary edge in normal direction;
% boundary_edges(1,k)=-2: Stress boundary edge in normal direction;
% boundary_edges(1,k)=-3: Robin boundary edge in normal direction.
% boundary_edges(2,k): specifiy the type of the kth boundary edge in tangential direction.
% boundary_edges(2,k)=-1: Dirichlet boundary edge in tangential direction;
% boundary_edges(2,k)=-2: Stress boundary edge in tangential direction;
% boundary_edges(2,k)=-3: Robin boundary edge in tangential direction.
% boundary_edges(3,k): index of the element which contains the kth boundary edge.
% boundary_edges(4:5,k): indices of the two end points of the kth boundary edge among all grid points, not the nodes of FE.
%                       That is, the index of partition is used here.
% boundary_edges(6:7,k): the unit outer normal vector at the kth boundary edge.
% boundary_edges(8:9,k): the unit tangential vector at the kth boundary edge. 
%                       It's counterclockwise to go from the normal vector to the tangential vector.

%% 2D Lagrange linear FE && 2D Lagrange quadratic FE
if type == 1 || type == 2
    xMeshStep = steps.x;
    yMeshStep = steps.y;
    
    if type == 2
        xEleStep = 2 * xMeshStep;
        yEleStep = 2 * yMeshStep;
    else
        xEleStep = xMeshStep;
        yEleStep = yMeshStep;
    end
    
    % boundary nodes
    obj.numNodes = 2 * (xEleStep + yEleStep);
    obj.nodes = zeros(obj.numNodes,7);
    
    obj.nodes(:,1) = -1;
    obj.nodes(:,2) = -1;
    
    % left node
    for i = 1 : yEleStep
        obj.nodes(i,3) = i;
        obj.nodes(i,4:5) = [-1,0];
    end
    
    % top node
    for i = yEleStep + 1 : (xEleStep + yEleStep)
        obj.nodes(i,3) = (i - yEleStep) * (yEleStep + 1);
        obj.nodes(i,4:5) = [0,1];
    end
    
    % right node
    for i = (xEleStep + yEleStep + 1) : (xEleStep + 2*yEleStep)
        obj.nodes(i,3) = (yEleStep + 1) * xEleStep + (i - xEleStep - yEleStep + 1);
        obj.nodes(i,4:5) = [1,0];
    end
    
    % bottom node
    for i = (xEleStep + 2*yEleStep + 1) : 2*(xEleStep + yEleStep)
        obj.nodes(i,3) = (i - xEleStep - 2*yEleStep) * (yEleStep + 1) + 1;
        obj.nodes(i,4:5) = [0,-1];
    end
    
    % Correct the normal direction at the four corners.
    % left-bottom
    obj.nodes(1,4:5) = [-1,-1]/sqrt(2);
    % left-top
    obj.nodes(yEleStep + 1,4:5) = [-1,1]/sqrt(2);
    % right-top
    obj.nodes(xEleStep + 2*yEleStep,4:5) = [1,1]/sqrt(2);
    % right-bottom
    obj.nodes(2*(xEleStep + yEleStep),4:5) = [1,-1]/sqrt(2);
    
    %It's counterclockwise to go from the normal vector n to the tangential vector \tau.
    %Hence \tau_1=-n_2, \tau_2=n_1.
    obj.nodes(:,6) = -obj.nodes(:,5);
    obj.nodes(:,7) = obj.nodes(:,4);
    
    % boundary edges
    obj.numEdges = 2 * (xMeshStep + yMeshStep);
    obj.edges = zeros(obj.numEdges,9);
    
    obj.edges(:,1) = -1;
    obj.edges(:,2) = -1;
    
    % left edge
    for i = 1 : yMeshStep
        obj.edges(i,3) = 2*i - 1;
        obj.edges(i,4) = i;
        obj.edges(i,5) = i + 1;
        obj.edges(i,6:7) = [-1,0];
    end
    
    % top edge
    for i = yMeshStep + 1 : xMeshStep + yMeshStep
        obj.edges(i,3) = (i - yMeshStep) * (2*yMeshStep);
        obj.edges(i,4) = (i - yMeshStep) * (yMeshStep + 1);
        obj.edges(i,5) = (i - yMeshStep + 1) * (yMeshStep + 1);
        obj.edges(i,6:7) = [0,1];
    end
    
    % right edge
    for i = xMeshStep + yMeshStep + 1 : xMeshStep + 2*yMeshStep
        obj.edges(i,3) = (xMeshStep-1)*(2*yMeshStep) + 2*(i-xMeshStep-yMeshStep);
        obj.edges(i,4) = xMeshStep * (yMeshStep + 1) + (i-xMeshStep-yMeshStep);
        obj.edges(i,5) = xMeshStep * (yMeshStep + 1) + (i-xMeshStep-yMeshStep) + 1;
        obj.edges(i,6:7) = [1,0];
    end
    
    % bottom edge
    for i = xMeshStep + 2*yMeshStep + 1 : 2*(xMeshStep + yMeshStep)
        obj.edges(i,3) = (i-xMeshStep-2*yMeshStep-1) * (2*yMeshStep) + 1;
        obj.edges(i,4) = (i-xMeshStep-2*yMeshStep-1) * (yMeshStep + 1) + 1;
        obj.edges(i,5) = (i-xMeshStep-2*yMeshStep) * (yMeshStep + 1) + 1;
        obj.edges(i,6:7) = [0,-1];
    end
    
    %It's counterclockwise to go from the normal vector n to the tangential vector \tau.
    %Hence \tau_1=-n_2, \tau_2=n_1.
    obj.edges(:,8) = -obj.edges(:,7);
    obj.edges(:,9) = obj.edges(:,6);
    
    obj.partition = [xMeshStep,yMeshStep,xEleStep,yEleStep];
    return
end

%% Warning
warning('Wrong TYPE of mesh!');
end

