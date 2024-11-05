classdef ElasticityTriSolver < handle
    
    properties
        meshInfo
        
        trialElementInfo
        testElementInfo
        
        numElements
        
        elementComputeMethod
        
        lhsMatrix
        rhsVector
        
        femSolution
        
        % TODO: boundaryProcess
        
        % other functions required to solve the problem
        fLambda
        fMu
        f1
        f2
    end
    
    methods
        function obj = ElasticityTriSolver(mesh,trialElement,testElement,trial,test,gaussQuad,lam,mu,g1,g2)
            obj.meshInfo = mesh;
            obj.trialElementInfo = trialElement;
            obj.testElementInfo = testElement;
            obj.numElements = mesh.numElements;
            obj.elementComputeMethod = TriElementCompute(trial,test,gaussQuad);
            obj.fLambda = lam;
            obj.fMu = mu;
            obj.f1 = g1;
            obj.f2 = g2;
        end
        
        generateEquation(obj);
        
        function solve(obj)
            obj.femSolution = obj.lhsMatrix\obj.rhsVector;
        end
        
    end
end

