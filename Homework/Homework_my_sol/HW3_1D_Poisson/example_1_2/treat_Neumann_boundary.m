function [A,b] = treat_Neumann_boundary(boundarynodes,boundaryvalues,coe_fun,N,A,b,mesh_type)

if mesh_type == 101
    A(1,:) = 0;
    A(1,1) = 1;
%     A(N,:) = 0;
%     A(N,N) = 1;
    b(1) = boundaryvalues(1);
    b(N) = b(N) + boundaryvalues(2)*coe_fun(boundarynodes(2));

end

end