classdef TriQuadratureTest < matlab.unittest.TestCase
   
    methods(Test)
        
        function simpleQuad(testCase)
            t_x = 1.3;
            t_y = 0.8;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            ymax = @(x) t_y*(1 - x/t_x);
            quadFun = @(x, y) x + 2*y;
            
            m_quad = TriQuadrature(9);
            actSolution = m_quad.evaluate(quadFun, vertices);
            expSolution = integral2(quadFun, 0, t_x, 0, ymax);
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-10);
        end
        
        function squareQuad(testCase)
            t_x = 0.8;
            t_y = 1.3;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            ymax = @(x) t_y.*(1 - x/t_x);
            quadFun = @(x, y) x.*x + 3*y.*y;
            
            m_quad = TriQuadrature(9);
            actSolution = m_quad.evaluate(quadFun, vertices);
            expSolution = integral2(quadFun, 0, t_x, 0, ymax);
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-10);
        end
        
        function exponentialQuad(testCase)
            % 没想到精度这么低
            t_x = 0.88;
            t_y = 1.88;
            vertices = [0, 0;
                        t_x, 0;
                        0, t_y];
            ymax = @(x) t_y.*(1 - x/t_x);
            quadFun = @(x, y) 2*exp(x) - exp(y);
            
            m_quad = TriQuadrature(9);
            actSolution = m_quad.evaluate(quadFun, vertices);
            expSolution = integral2(quadFun, 0, t_x, 0, ymax);
            testCase.verifyEqual(actSolution,expSolution,"AbsTol",1e-4);
        end
        
    end
end

