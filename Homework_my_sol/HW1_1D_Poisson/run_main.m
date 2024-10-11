clear
clc

%% 
x_domain = [0,1];
% 求解区域

boundaryvalues = [0,cos(1)];
% 边界值

N = 4;
% N = 8;
% 划分单元个数

mesh_type = 101;
% 单元类型
% 101: 1D linear

basis_type_trial = 101;
% basis_type_trial = 102;
% trial function 
% 101: 1D linear

basis_type_test = 101;
% basis_type_test = 102;
% test function
% 101: 1D linear

% coe
% coefficient: c(x)

% load
% load: f(x)

Gauss_type = 104;
% 高斯求积类型
% Gauss_type == 10i: 1表示1维，i表示求积点个数

% ref
% 参考解

%% solver
[solution,error] = FE_solver_1D_Poisson(x_domain,boundaryvalues,N,mesh_type,basis_type_trial,basis_type_test,@coe,@load,@ref,Gauss_type);

%% plot
fprintf('%.4e\n',error);
%% function
function coefficient = coe(x)

coefficient = exp(x);

end

function load_fun = load(x)

load_fun = -exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x));

end

function reference = ref(x)

reference = x*cos(x);

end



















