classdef SteadyStokesTriSolver < handle

    properties
        % Geometry
        meshInfo
        
        uTrialElementInfo
        uTestElementInfo
        pTrialElementInfo
        pTestElementInfo
        
        numElements
        
        boundaryInfo
        
        % Equation
        uuComputeMethod
        puComputeMethod
        
        lhsMatrix
        rhsVector
        
        femSolution
        
        % other functions required to solve the problem
        fMu
        f1
        f2
        boundF1
        boundF2
    end
    
    methods
        function obj = SteadyStokesTriSolver(mesh,uTrialElement,uTestElement,pTrialElement,pTestElement,boundary,uTrial,uTest,pTrial,pTest,gauss,mu,g1,g2,bf1,bf2)
            obj.meshInfo = mesh;
            obj.uTrialElementInfo = uTrialElement;
            obj.uTestElementInfo = uTestElement;
            obj.pTrialElementInfo = pTrialElement;
            obj.pTestElementInfo = pTestElement;
            obj.numElements = mesh.numElements;
            obj.boundaryInfo = boundary;
            
            obj.uuComputeMethod = TriElementCompute(uTrial,uTest,gauss);
            obj.puComputeMethod = TriElementCompute(pTrial,uTest,gauss);
            
            obj.fMu = mu;
            obj.f1 = g1;
            obj.f2 = g2;
            obj.boundF1 = bf1;
            obj.boundF2 = bf2;
        end
        
        generateEquations(obj);
        
        boundaryConditions(obj);
        
        function solve(obj)
            uColSize = obj.uTrialElementInfo.numPoints;
            pColSize = obj.pTrialElementInfo.numPoints;
            x = obj.lhsMatrix\obj.rhsVector;
            solution.u1 = x(1:uColSize);
            solution.u2 = x(uColSize+1:2*uColSize);
            solution.p = x(2*uColSize+1:2*uColSize+pColSize);
            obj.femSolution = solution;
        end
    end
end

