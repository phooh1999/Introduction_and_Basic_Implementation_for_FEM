function [A,b] = treat_Robin_boundary(boundarynodes,boundaryvalues,coe_fun,N,A,b,mesh_type)

if mesh_type == 101
%     A(1,:) = 0;
%     A(1,1) = 1;
    A(1,1) = A(1,1) - coe_fun(boundarynodes(1));
    A(N,:) = 0;
    A(N,N) = 1;
    b(1) = b(1)-coe_fun(boundarynodes(1))*boundaryvalues(1);
    b(N) = boundaryvalues(2);
% 这里写得比较case by case，等二维再系统处理边界条件
end

end