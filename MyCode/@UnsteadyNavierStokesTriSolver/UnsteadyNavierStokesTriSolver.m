classdef UnsteadyNavierStokesTriSolver < handle

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
        
        MMatrix
        AMatrix
        
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
        function obj = UnsteadyNavierStokesTriSolver(mesh,uTrialElement,uTestElement,pTrialElement,pTestElement,boundary,uTrial,uTest,pTrial,pTest,gauss,mu,g1,g2,bf1,bf2)
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
        
        generateMatrix(obj);
        
        bVec = generateVector(obj,t);
        
        [AN, bN] = generateNonlinearPart(obj,preSol);
        
        boundaryConditions(obj,t);
        
        fixPressure(obj,funRefp,index,t);
        
        x0 = generateInitialVal(obj,fu1,fu2,fp);
        
        solve(obj,x0,t0,t1,dt,theta,maxStep,rtol,atol);

    end
end

