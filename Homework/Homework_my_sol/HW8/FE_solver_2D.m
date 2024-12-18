function [solution] = FE_solver_2D(x_domain,lambda_fun,mu_fun,bound_fun,N,mesh_type,basis_type_trial,basis_type_test,load_fun_1,load_fun_2,Gauss_type)
% mesh_type == 201: 2D uniform partition
% basis_type_~ == 201: 2D linear
% basis_type_~ == 202: 2D quadratic

[P,T] = generate_PT(x_domain,N,mesh_type);
boundaryedges = generate_boundaryedges(N,T);
% 划分单元网格，获取边界边
% x_domain: 求解区域
% N划分网格数

if basis_type_trial == 201
    Pb_trial = P;
    Tb_trial = T;
    number_of_local_basis_fun_trial = 3;
elseif basis_type_trial == 202
    [Pb_trial,Tb_trial] = generate_PT(x_domain,N,202);
    number_of_local_basis_fun_trial = 6;
end

if basis_type_test == 201
    Pb_test = P;
    Tb_test = T;
    number_of_local_basis_fun_test = 3;
elseif basis_type_test == 202
    [Pb_test,Tb_test] = generate_PT(x_domain,N,202);
    number_of_local_basis_fun_test = 6;
end

[boundaryedges1] = generate_boundaryedges(N,Tb_trial);
[boundarynodes] = generate_boundarynodes(boundaryedges1,basis_type_trial);
boundarynodes(1,1) = -1;
% 先不固定，试一下
% 确实没有固定成dirichlet的准
clear boundaryedges1
% 获取边界节点
% 由于之前理解不到位，没有将网格和有限元单元区分开来
% 这里强行指定为Dirichlet条件吧，看之后参考程序怎么写的

matrix_size = [size(Pb_test,2),size(Pb_trial,2)];
% 列为trial的未知数个数，目前为点个数
% 行为test的个数，目前也为点个数

number_of_elements = size(T,2);
% T的第2维就是单元个数

A1 = assemble_matrix_2D(lambda_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,1,0,basis_type_test,1,0);
A2 = assemble_matrix_2D(mu_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,1,0,basis_type_test,1,0);
A3 = assemble_matrix_2D(mu_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,0,1,basis_type_test,0,1);
A4 = assemble_matrix_2D(lambda_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,0,1,basis_type_test,1,0);
A5 = assemble_matrix_2D(mu_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,1,0,basis_type_test,0,1);
A6 = assemble_matrix_2D(lambda_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,1,0,basis_type_test,0,1);
A7 = assemble_matrix_2D(mu_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,0,1,basis_type_test,1,0);
A8 = assemble_matrix_2D(lambda_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,0,1,basis_type_test,0,1);

A = [A1+2*A2+A3,A4+A5;A6+A7,A8+2*A3+A2];

b1 = assemble_vector_2D(load_fun_1,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);
b2 = assemble_vector_2D(load_fun_2,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);

b = [b1;b2];

[A,b] = treat_Dirichlet_boundary_stress(boundarynodes,bound_fun,Pb_trial,A,b);

solution = A\b;

%% error

% error = FE_max_err_compute(ref_fun,Pb_trial,solution);
% [error] = compute_Hs_error(solution,ref_fun,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,der_s);
% [error] = compute_L_infi_error(solution,ref_fun,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,der_s);

end
