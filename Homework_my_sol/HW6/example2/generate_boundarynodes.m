function [boundarynodes] = generate_boundarynodes(boundaryedges,mesh_type)

% GENERATE_BOUNDARYNODES 此处显示有关此函数的摘要
% 此处显示详细说明

if mesh_type == 201
    
    nbn = size(boundaryedges,2);
    boundarynodes = zeros(2,nbn);
    boundarynodes(1,:) = boundaryedges(1,:);
    boundarynodes(2,:) = boundaryedges(3,:);
    
elseif mesh_type == 202
    
    nbn = size(boundaryedges,2);
    boundarynodes = zeros(2,2*nbn);
    
    for i = 1:nbn
        
        boundarynodes(1,2*i-1) = boundaryedges(1,i);
        boundarynodes(2,2*i-1) = boundaryedges(3,i);
        boundarynodes(1,2*i) = boundaryedges(1,i);
        boundarynodes(2,2*i) = (boundaryedges(3,i)+boundaryedges(4,i))/2;
        
    end
    
end

end

