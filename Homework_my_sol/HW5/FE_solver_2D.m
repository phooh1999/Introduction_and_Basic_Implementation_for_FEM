function [solution] = FE_solver_2D(x_domain,bound_fun,N,mesh_type,basis_type_trial,basis_type_test,coe_fun,load_fun,Gauss_type)
% mesh_type == 201: 2D uniform partition
% basis_type_~ == 201: 2D linear
% basis_type_~ == 202: 2D quadratic

[P,T] = generate_PT(x_domain,N,mesh_type);
boundaryedges = generate_boundaryedges(N,T);
% ���ֵ�Ԫ���񣬻�ȡ�߽��
% x_domain: �������
% N����������

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
clear boundaryedges1
% ��ȡ�߽�ڵ�
% ����֮ǰ���ⲻ��λ��û�н����������Ԫ��Ԫ���ֿ���

matrix_size = [size(Pb_test,2),size(Pb_trial,2)];
% ��Ϊtrial��δ֪��������ĿǰΪ�����
% ��Ϊtest�ĸ�����ĿǰҲΪ�����

number_of_elements = size(T,2);
% T�ĵ�2ά���ǵ�Ԫ����

basis_der_x_trial = 1;
basis_der_x_test = 1;
basis_der_y_trial = 0;
basis_der_y_test = 0;
A1 = assemble_matrix_2D(coe_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_der_x_test,basis_der_y_test);

basis_der_x_trial = 0;
basis_der_x_test = 0;
basis_der_y_trial = 1;
basis_der_y_test = 1;
A2 = assemble_matrix_2D(coe_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_der_x_test,basis_der_y_test);

A = A1 + A2;
% 2D����Ҫ����������װ����

load_basis_der_x_test = 0;
load_basis_der_y_test = 0;
% ע��load����Ĳ��Ժ���û����
b = assemble_vector_2D(load_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,load_basis_der_x_test,load_basis_der_y_test);

[A,b] = treat_boundary(boundarynodes,bound_fun,Pb_trial,A,b);

solution = A\b;

%% error

% error = FE_max_err_compute(ref_fun,Pb_trial,solution);
% [error] = compute_Hs_error(solution,ref_fun,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,der_s);
% [error] = compute_L_infi_error(solution,ref_fun,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,der_s);

end