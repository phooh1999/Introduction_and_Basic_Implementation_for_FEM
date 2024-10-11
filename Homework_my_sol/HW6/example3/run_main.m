clear
% clc

%% 
x_domain = [-1,1,-1,1];
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

for n = [16 32 64 128]
%% solver
N = [n,n];
[solution] = FE_solver_2D(x_domain,@fun_robin_r,@fun_robin_q,@fun_neumann,@bound,N,mesh_type,basis_type_trial,basis_type_test,@coe,@load,Gauss_type);
% 需要把计算误差的函数拿到外面来，减少重复计算
[P,T] = generate_PT(x_domain,N,mesh_type);
[Pb_trial,Tb_trial] = generate_PT(x_domain,N,basis_type_trial);
if basis_type_trial == 201
    number_of_local_basis_fun_trial = 3;
elseif basis_type_trial == 202
    number_of_local_basis_fun_trial = 6;
end
% error_max = FE_max_err_compute(@ref,Pb_trial,solution);
error_infi = compute_L_infi_error(solution,@ref,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_h0 = compute_Hs_error(solution,@ref,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,0);
error_h1 = compute_Hs_error(solution,@ref,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun_trial,basis_type_trial,1);
%% print
fprintf('N = %d\n',n);
% fprintf('%.4e\n',error_max);
fprintf('%.4e\n',error_infi);
fprintf('%.4e\n',error_h0);
fprintf('%.4e\n',error_h1);
fprintf('\n');
end

%% function
function coefficient = coe(x,y)

coefficient = 1;

end

function result = fun_neumann(x,y)

result = -exp(x-1);

end

function r = fun_robin_r(x,y)

r = 1;

end

function q = fun_robin_q(x,y)

q = 0;

end

function coefficient = bound(x,y)

if y == 1
    coefficient = exp(x+1);
%     coefficient = 0;
% elseif y == -1
%     coefficient = -exp(x-1);
% %     coefficient = -2*x*(1-x/2)*exp(x-1);
elseif x == 1
    coefficient = exp(1+y);
%     coefficient = 0.5*y*(1-y)*exp(1+y);
elseif x == -1
    coefficient = exp(-1+y);
%     coefficient = -1.5*y*(1-y)*exp(-1+y);
    
end

end

function load_fun = load(x,y)

% load_fun = -y*(1-y)*(1-x-x^2/2)*exp(x+y)-x*(1-x/2)*(-3*y-y^2)*exp(x+y);
load_fun = -2*exp(x+y);

end

function reference = ref(x,y,basis_der_x,basis_der_y)

if basis_der_x == 0 && basis_der_y == 0
    reference = exp(x+y);
%     reference = x*y*(1-x/2)*(1-y)*exp(x+y);
elseif basis_der_x == 1 && basis_der_y == 0
    reference = exp(x+y);
%     reference = y*(1-y)*(1-x^2/2)*exp(x+y);
elseif basis_der_x == 0 && basis_der_y == 1
    reference = exp(x+y);
%     reference = x*(1-x/2)*(1-y-y^2)*exp(x+y);
end

end