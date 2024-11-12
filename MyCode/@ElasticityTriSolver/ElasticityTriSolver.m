classdef ElasticityTriSolver < handle
    
    properties
        % Geometry
        meshInfo
        
        trialElementInfo
        testElementInfo
        
        numElements
        
        boundaryInfo
        
        % Equation
        elementComputeMethod
        
        lhsMatrix
        rhsVector
        
        femSolution
        
        % other functions required to solve the problem
        fLambda
        fMu
        f1
        f2
        boundF1
        boundF2
    end
    
    methods
        function obj = ElasticityTriSolver(mesh,trialElement,testElement,boundary,trial,test,gauss,lam,mu,g1,g2,bf1,bf2)
            obj.meshInfo = mesh;
            obj.trialElementInfo = trialElement;
            obj.testElementInfo = testElement;
            obj.numElements = mesh.numElements;
            obj.boundaryInfo = boundary;
            obj.elementComputeMethod = TriElementCompute(trial,test,gauss);
            obj.fLambda = lam;
            obj.fMu = mu;
            obj.f1 = g1;
            obj.f2 = g2;
            obj.boundF1 = bf1;
            obj.boundF2 = bf2;
        end
        
        generateEquations(obj);
        
        boundaryConditions(obj);
        
        function solve(obj)
            colSize = obj.trialElementInfo.numPoints;
            obj.femSolution = zeros(colSize,2);
            x = obj.lhsMatrix\obj.rhsVector;
            solution.u1 = x(1:colSize);
            solution.u2 = x(colSize+1:2*colSize);
            obj.femSolution = solution;
        end
        
    end
end

