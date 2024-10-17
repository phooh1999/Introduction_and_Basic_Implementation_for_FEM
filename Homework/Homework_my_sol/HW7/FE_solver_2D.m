function [solution] = FE_solver_2D(x_domain,dt,t_span,theta,fun_ini,fun_const_1,bound_fun,N,mesh_type,basis_type_trial,basis_type_test,coe_fun,load_fun,Gauss_type)
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
clear boundaryedges1
% 获取边界节点
% 由于之前理解不到位，没有将网格和有限元单元区分开来

matrix_size = [size(Pb_test,2),size(Pb_trial,2)];
% 列为trial的未知数个数，目前为点个数
% 行为test的个数，目前也为点个数

number_of_elements = size(T,2);
% T的第2维就是单元个数

M = assemble_matrix_2D(fun_const_1,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,0,0,basis_type_test,0,0);

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
% 2D问题要调用两次组装函数

A_tilde = M/dt+theta*A;
A_fixed = M/dt-(1-theta)*A;
clear M A A1 A2

t_step = (t_span(1,2)-t_span(1,1))/dt;

u_m = generate_initial_value(Pb_trial,fun_ini);

load_basis_der_x_test = 0;
load_basis_der_y_test = 0;
% 注意load这里的测试函数没有求导

for i = 0:(t_step-1)
    
    t_mp1 = (i+1)*dt;
    t_m = i*dt;
    b_mp1 = assemble_vector_2D_t(t_mp1,load_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,load_basis_der_x_test,load_basis_der_y_test);
    b_m = assemble_vector_2D_t(t_m,load_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,load_basis_der_x_test,load_basis_der_y_test);
    
    b_tilde = theta*b_mp1+(1-theta)*b_m+A_fixed*u_m;
    
    [A,b] = treat_Dirichlet_boundary_t(t_mp1,boundarynodes,bound_fun,Pb_trial,A_tilde,b_tilde);
    
    u_mp1 = A\b;
    
    solution = u_mp1;
    
    u_m = u_mp1;
    
end

end
