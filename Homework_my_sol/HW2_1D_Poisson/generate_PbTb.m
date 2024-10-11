function [Pb,Tb] = generate_PbTb(x,N,mesh_type)
% x: 划分区域
% N: 划分网格数
% mesh_type: 网格划分类型
% mesh_type == 101: 1D uniform partition
%% 一维网格
if mesh_type == 101
    Pb = zeros(1,N+1);
    h = (x(2) - x(1))/N;
    Pb(1,:) = x(1):h:x(2);
    
    Tb = zeros(2,N);
    Tb(1,:) = 1:N;
    Tb(2,:) = 2:(N+1);
%% 单元加中点
elseif mesh_type == 102
    Pb = zeros(1,2*N+1);
    h = (x(2) - x(1))/(2*N);
    Pb(1,:) = x(1):h:x(2);
    
    Tb = zeros(3,N);
    Tb(1,:) = 1:2:(2*N-1);
    Tb(2,:) = 3:2:(2*N+1);
    Tb(3,:) = 2:2:2*N;
%   Tb(1,k)第k个单元的左端点
%   Tb(2,k)第k个单元的右端点
%   Tb(3,k)第k个单元的中点
end

end