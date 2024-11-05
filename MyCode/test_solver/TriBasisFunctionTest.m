classdef TriBasisFunctionTest < matlab.unittest.TestCase
   
    methods(Test)
        
        function constant2D(testCase)
            t_x = 1.3;
            t_y = 0.8;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            
            basisFun = LocalTriBasisFunction(0);
            
            actSolution = basisFun.evaluate(0.2,0.3,vertices,0,0);
            expSolution = 1;
            testCase.verifyEqual(actSolution,expSolution);
            
        end
        
        function lagrangeLinear2D(testCase)
            t_x = 1.3;
            t_y = 0.8;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            
            basisFun = LocalTriBasisFunction(1);
            
            actSolution = basisFun.evaluate(0.2,0.3,vertices,0,0);
            expSolution = [0.471153846153846;
                           0.153846153846154;
                           0.375000000000000];
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-15);
            
        end
        
        function lagrangeQuadratic2D(testCase)
            t_x = 1.3;
            t_y = 0.8;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            
            basisFun = LocalTriBasisFunction(2);
            
            actSolution = basisFun.evaluate(0.2,0.3,vertices,0,0);
            expSolution = [-0.027181952662722;
                           -0.106508875739645;
                           -0.093750000000000;
                            0.289940828402367;
                            0.230769230769231;
                            0.706730769230769;
                           ];
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-15);
            
            actSolution = basisFun.evaluate(0.2,0.3,vertices,1,0);
            expSolution = [-0.680473372781065;
                           -0.295857988165680;
                            0;
                            0.976331360946745;
                            1.153846153846154;
                           -1.153846153846154
                           ];
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-15);
            
            actSolution = basisFun.evaluate(0.2,0.3,vertices,0,1);
            expSolution = [-1.105769230769231;
                            0;
                            0.625000000000000;
                           -0.769230769230770;
                            0.769230769230770;
                            0.480769230769231
                           ];
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-15);
            
        end
        
    end
end

