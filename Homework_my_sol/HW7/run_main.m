clear
% clc

%% 
x_domain = [0,2,0,1];
% 求解区域

t_span = [0,1];

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

theta = 0.5;

% for n = [8 16 32 64 128]
%% solver
n = 32;
N = [n,n/2];
dt = 1/64;
[solution] = FE_solver_2D(x_domain,dt,t_span,theta,@inni,@fun_1,@bound,N,mesh_type,basis_type_trial,basis_type_test,@coe,@load,Gauss_type);
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
% end

%% function
function coefficient = coe(x,y,t)

coefficient = 2;

end

function coefficient = fun_1(x,y)

coefficient = 1;

end

function initial_value = inni(x,y)

initial_value = exp(x+y);

end

function coefficient = bound(x,y,t)

if x == 0
    coefficient = exp(y+t);
%     coefficient = 0;
elseif x == 2
    coefficient = exp(2+y+t);
% %     coefficient = -2*x*(1-x/2)*exp(x-1);
elseif y == 0
    coefficient = exp(x+t);
%     coefficient = 0.5*y*(1-y)*exp(1+y);
elseif y == 1
    coefficient = exp(x+1+t);
%     coefficient = -1.5*y*(1-y)*exp(-1+y);
    
end

end

function load_fun = load(x,y,t)

load_fun = -3*exp(x+y+t);

end

function reference = ref(x,y,basis_der_x,basis_der_y)

if basis_der_x == 0 && basis_der_y == 0
    reference = exp(x+y+1);
%     reference = x*y*(1-x/2)*(1-y)*exp(x+y);
elseif basis_der_x == 1 && basis_der_y == 0
    reference = exp(x+y+1);
%     reference = y*(1-y)*(1-x^2/2)*exp(x+y);
elseif basis_der_x == 0 && basis_der_y == 1
    reference = exp(x+y+1);
%     reference = x*(1-x/2)*(1-y-y^2)*exp(x+y);
end

end