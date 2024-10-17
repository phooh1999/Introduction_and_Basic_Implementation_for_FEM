function [A,b] = treat_Dirichlet_boundary(boundarynodes,boundaryvalues,N,A,b,mesh_type)

if mesh_type == 101
    A(1,:) = 0;
    A(1,1) = 1;
    A(N+1,:) = 0;
    A(N+1,N+1) = 1;
    b(1) = boundaryvalues(1);
    b(N+1) = boundaryvalues(2);
end
% 一维问题先这么写，之后修改

end