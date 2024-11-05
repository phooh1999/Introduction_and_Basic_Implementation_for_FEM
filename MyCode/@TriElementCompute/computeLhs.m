function A = computeLhs(obj,vertices,coe,trialOrderX,trialOrderY,testOrderX,testOrderY)

% attention!!!
% here we use vector
% because A_ij = trial_j * test_i
% we should compute test * trial'
fun = @(x,y)  coe(x,y)...
            * obj.testFun.evaluate(x,y,vertices,testOrderX,testOrderY)...
            * (obj.trialFun.evaluate(x,y,vertices,trialOrderX,trialOrderY))';
        
A = obj.quadrature2D.evaluate(fun,vertices);

end