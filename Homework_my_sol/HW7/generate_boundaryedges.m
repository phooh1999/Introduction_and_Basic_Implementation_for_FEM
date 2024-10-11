function [boundaryedges] = generate_boundaryedges(N,T)

% GENERATE_BOUNDARYEDGES 
% -1 Dirichlet
% -2 Neumann

%% mesh_type = 201;

number_of_edges = 2*(N(1)+N(2));
boundaryedges = zeros(4,number_of_edges);

% Bottom

boundaryedges(1,1:N(1)) = -1;
% boundaryedges(1,1:N(1)) = -2;
% boundaryedges(1,1:N(1)) = -3;
boundaryedges(2,1:N(1)) = (0:(N(1)-1))*2*N(2)+1;
boundaryedges(3,1:N(1)) = T(1,boundaryedges(2,1:N(1)));
boundaryedges(4,1:N(1)) = T(2,boundaryedges(2,1:N(1)));

% right

boundaryedges(1,N(1)+1:N(1)+N(2)) = -1;
boundaryedges(2,N(1)+1:N(1)+N(2)) = boundaryedges(2,N(1))+1:2:boundaryedges(2,N(1))+1+2*(N(2)-1);
boundaryedges(3,N(1)+1:N(1)+N(2)) = T(2,boundaryedges(2,N(1)+1:N(1)+N(2)));
boundaryedges(4,N(1)+1:N(1)+N(2)) = T(3,boundaryedges(2,N(1)+1:N(1)+N(2)));

% top

boundaryedges(1,N(1)+N(2)+1:2*N(1)+N(2)) = -1;
boundaryedges(2,N(1)+N(2)+1:2*N(1)+N(2)) = boundaryedges(2,N(1)+N(2)):-2*N(2):boundaryedges(2,N(1)+N(2))-2*N(2)*(N(1)-1);
boundaryedges(3,N(1)+N(2)+1:2*N(1)+N(2)) = T(3,boundaryedges(2,N(1)+N(2)+1:2*N(1)+N(2)));
boundaryedges(4,N(1)+N(2)+1:2*N(1)+N(2)) = T(1,boundaryedges(2,N(1)+N(2)+1:2*N(1)+N(2)));

% left

boundaryedges(1,2*N(1)+N(2)+1:2*(N(1)+N(2))) = -1;
boundaryedges(2,2*N(1)+N(2)+1:2*(N(1)+N(2))) = boundaryedges(2,2*N(1)+N(2))-1:-2:boundaryedges(2,1);
boundaryedges(3,2*N(1)+N(2)+1:2*(N(1)+N(2))) = T(3,boundaryedges(2,2*N(1)+N(2)+1:2*(N(1)+N(2))));
boundaryedges(4,2*N(1)+N(2)+1:2*(N(1)+N(2))) = T(1,boundaryedges(2,2*N(1)+N(2)+1:2*(N(1)+N(2))));

end

