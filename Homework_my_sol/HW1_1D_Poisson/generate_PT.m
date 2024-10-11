function [P,T] = generate_PT(x,N,mesh_type)
% x: 划分区域
% N: 划分网格数
% mesh_type: 网格划分类型
% mesh_type == 101: 1D uniform partition

if mesh_type == 101
    P = zeros(1,N+1);
    h = (x(2) - x(1))/N;
    P(1,:) = x(1):h:x(2);
    
    T = zeros(2,N);
    T(1,:) = 1:N;
    T(2,:) = 2:(N+1);
end

end