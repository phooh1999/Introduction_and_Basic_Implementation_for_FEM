%% ---- 测试... ----

%% 测试三角形高斯积分
disp('----测试三角形高斯积分----')
fprintf('\n')

testCase1 = TriQuadratureTest;
testCase1.run

fprintf('\n')

%% 测试局部基函数
disp('----测试局部基函数----')
fprintf('\n')

testCase2 = TriBasisFunctionTest;
testCase2.run

fprintf('\n')