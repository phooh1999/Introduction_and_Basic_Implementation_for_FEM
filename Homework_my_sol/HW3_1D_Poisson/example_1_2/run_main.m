clear
clc

%% 
x_domain = [0,1];
% 求解区域

% boundaryvalues = [0,cos(1)]; %Drichlet
boundaryvalues = [0,cos(1)-sin(1)]; % Neumann
% boundaryvalues = [1,cos(1)]; %Robin

% 边界值

% N = 4;
% N = 8;
% 划分单元个数

mesh_type = 101;
% 网格单元类型
% 101: 1D linear

% basis_type_trial = 101;
basis_type_trial = 102;
% trial function 
% 101: 1D linear

% basis_type_test = 101;
basis_type_test = 102;
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

% error
% 误差计算 der_s = 0,1
% L2 norm(s=0), H1 semi-norm(s=1), L Infinity norm(s=0,L_infi_error)
% der_s = 0;
der_s = 1;

for N = [4 8 16 32 64 128]
%% solver
[solution,error] = FE_solver_1D_Poisson(x_domain,boundaryvalues,N,mesh_type,basis_type_trial,basis_type_test,@coe,@load,@ref,Gauss_type,der_s);

%% print
fprintf('%.4e\n',error);
end

%% function
function coefficient = coe(x)

coefficient = exp(x);

end

function load_fun = load(x)

load_fun = -exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x));

end

function reference = ref(x,basis_der_x)

if basis_der_x == 0
    reference = x*cos(x);
elseif basis_der_x == 1
    reference = cos(x)-x*sin(x);
end

end



















