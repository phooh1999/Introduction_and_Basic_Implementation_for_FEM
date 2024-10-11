function [solution,error] = FE_solver_1D_Poisson(x_domain,boundaryvalues,N,mesh_type,basis_type_trial,basis_type_test,coe_fun,load_fun,ref_fun,Gauss_type)
% mesh_type == 101: 1D uniform partition
% basis_type_~ == 101: 1D linear

[P,T] = generate_PT(x_domain,N,mesh_type);
% 一维问题P,T和Pb,Tb是一致的，暂时不做区分
% x_domain: 求解区域
% N划分网格数

if basis_type_trial == 101
    Pb_trial = P;
    Tb_trial = T;
    number_of_local_basis_fun_trial = 2;
% elseif basis_type_trial == 102
%     Pb_trial = ?;
%     Tb_trial = ?;
%     number_of_local_basis_fun_trial = 3;
%     暂时不写二阶单元
% 在这个问题中，有限元节点和单元节点是一样的
% 先这么写，后续再研究
end

if basis_type_test == 101
    Pb_test = P;
    Tb_test = T;
    number_of_local_basis_fun_test = 2;
% elseif basis_type_test == 102
%     ?;
% 二阶问题暂时不考虑
% 取一样的基函数？这里先这么写
end

boundarynodes = generate_boundarynodes(x_domain,mesh_type);
% 获取边界节点
% 目前只写一维，即计算区域两端

matrix_size = [size(Pb_test,2),size(Pb_trial,2)];
% 初步猜测
% 列为trial的未知数个数，目前为点个数
% 行为test的个数，目前也为点个数，后续需要再次考虑

number_of_elements = size(T,2);
% T的第2维就是单元个数

basis_der_x_trial = 1;
basis_der_x_test = 1;
load_basis_der_x_test = 0;
% 注意load这里的测试函数没有求导

A = assemble_matrix_1D(coe_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test);

b = assemble_vector_1D(load_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,load_basis_der_x_test);

[A,b] = treat_Dirichlet_boundary(boundarynodes,boundaryvalues,number_of_elements,A,b,mesh_type);

solution = A\b;

%% error

error = FE_max_err_compute(ref_fun,Pb_trial,solution);

end








