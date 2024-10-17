clear
% clc

%% 
x_domain = [0,1,0,1];
% 求解区域

mesh_type = 201;
% 网格单元类型

% basis_type_trial = 201;
basis_type_trial = 202;
% trial function 

% basis_type_test = 201;
basis_type_test = 202;
% test function

% coe
% coefficient: c(x)

% load
% load: f(x)

Gauss_type = 209;
% 高斯求积类型

% ref
% 参考解

% error
% 误差计算 der_s = 0,1
% L2 norm(s=0), H1 semi-norm(s=1), L Infinity norm(s=0,L_infi_error)
% der_s = 0;
% der_s = 1;

for n = [8 16 32 64]
%% solver
N = [n,n];
[solution] = FE_solver_2D(x_domain,@lambda_fun,@mu_fun,@bound,N,mesh_type,basis_type_trial,basis_type_test,@load_fun_1,@load_fun_2,Gauss_type);
% 需要把计算误差的函数拿到外面来，减少重复计算
[P,T] = generate_PT(x_domain,N,mesh_type);
[Pb_trial,Tb_trial] = generate_PT(x_domain,N,basis_type_trial);
if basis_type_trial == 201
    number_of_local_basis_fun_trial = 3;
elseif basis_type_trial == 202
    number_of_local_basis_fun_trial = 6;
end

number_of_solution_1 = size(solution,1)/2;
number_of_solution_2 = size(solution,1);
solution1(:,1) = solution(1:number_of_solution_1,1);
solution2(:,1) = solution(number_of_solution_1+1:number_of_solution_2,1);

error_infi_1 = compute_L_infi_error(solution1,@ref_1,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_infi_2 = compute_L_infi_error(solution2,@ref_2,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_infi = max(error_infi_1,error_infi_2);

error_h0_1 = compute_Hs_error(solution1,@ref_1,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_h0_2 = compute_Hs_error(solution2,@ref_2,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_h0 = sqrt(error_h0_1^2+error_h0_2^2);

error_h1_1 = compute_Hs_error(solution1,@ref_1,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,1);
error_h1_2 = compute_Hs_error(solution2,@ref_2,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,1);
error_h1 = sqrt(error_h1_1^2+error_h1_2^2);

clear solution1 solution2
%% print
fprintf('N = %d\n',n);
% fprintf('%.4e\n',error_max);
fprintf('%.4e\n',error_infi);
fprintf('%.4e\n',error_h0);
fprintf('%.4e\n',error_h1);
fprintf('\n');
end

%% function

function [coefficient1,coefficient2] = bound(x,y)

if y == 0
    coefficient1 = 0;
    coefficient2 = 0;
elseif y == 1
    coefficient1 = 0;
    coefficient2 = 0;
elseif x == 0
    coefficient1 = 0;
    coefficient2 = 0;
elseif x == 1
    coefficient1 = 0;
    coefficient2 = 0;
end

end

function lambda = lambda_fun(x,y)
lambda = 1;
end

function mu = mu_fun(x,y)
mu = 2;
end

function [load_1] = load_fun_1(x,y)
lambda = 1;
mu = 2;

load_1 = -(lambda+2*mu)*(-pi^2*sin(pi*x)*sin(pi*y))-(lambda+mu)*((2*x-1)*(2*y-1))-mu*(-pi^2*sin(pi*x)*sin(pi*y));
% load_2 = -(lambda+2*mu)*(2*x*(x-1))-(lambda+mu)*(pi^2*cos(pi*x)*cos(pi*y))-mu*(2*y*(y-1));

end

function [load_2] = load_fun_2(x,y)
lambda = 1;
mu = 2;

% load_1 = -(lambda+2*mu)*(-pi^2*sin(pi*x)*sin(pi*y))-(lambda+mu)*((2*x-1)*(2*y-1))-mu*(-pi^2*sin(pi*x)*sin(pi*y));
load_2 = -(lambda+2*mu)*(2*x*(x-1))-(lambda+mu)*(pi^2*cos(pi*x)*cos(pi*y))-mu*(2*y*(y-1));

end

function [reference1] = ref_1(x,y,basis_der_x,basis_der_y)

if basis_der_x == 0 && basis_der_y == 0
    reference1 = sin(pi*x)*sin(pi*y);
%     reference2 = x*(x-1)*y*(y-1);
elseif basis_der_x == 1 && basis_der_y == 0
    reference1 = pi*cos(pi*x)*sin(pi*y);
%     reference2 = (2*x-1)*y*(y-1);
elseif basis_der_x == 0 && basis_der_y == 1
    reference1 = pi*sin(pi*x)*cos(pi*y);
%     reference2 = x*(x-1)*(2*y-1);
end

end

function [reference2] = ref_2(x,y,basis_der_x,basis_der_y)

if basis_der_x == 0 && basis_der_y == 0
%     reference1 = sin(pi*x)*sin(pi*y);
    reference2 = x*(x-1)*y*(y-1);
elseif basis_der_x == 1 && basis_der_y == 0
%     reference1 = pi*cos(pi*x)*sin(pi*y);
    reference2 = (2*x-1)*y*(y-1);
elseif basis_der_x == 0 && basis_der_y == 1
%     reference1 = pi*sin(pi*x)*cos(pi*y);
    reference2 = x*(x-1)*(2*y-1);
end

end