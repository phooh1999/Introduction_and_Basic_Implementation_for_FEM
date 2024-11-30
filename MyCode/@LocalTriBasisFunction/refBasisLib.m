function refBasisLib(obj,type)

if type == 0
    obj.numBasis = 1;
    obj.refBasisFun = @constant2D;
end

if type == 1
    obj.numBasis = 3;
    obj.refBasisFun = @lagrangeLinear2D;
end

if type == 2
    obj.numBasis = 6;
    obj.refBasisFun = @lagrangeQuadratic2D;
end

if type == 3
    obj.numBasis = 3;
    obj.refBasisFun = @CrouzeixRaviart2D;
end

end

function f = constant2D(~,~,orderX,orderY)

if orderX == 0 && orderY == 0
    f = 1;
end

end

function f = lagrangeLinear2D(x,y,orderX,orderY)

if orderX == 0 && orderY == 0
    
    f = [1-x-y; x; y];

elseif orderX == 1 && orderY == 0
    
    f = [-1; 1; 0];

elseif orderX == 0 && orderY == 1
    
    f = [-1; 0; 1];
        
end

end

function f = lagrangeQuadratic2D(x,y,orderX,orderY)

if orderX == 0 && orderY == 0

    f = [1-3*x-3*y+2*x.^2+2*y.^2+4*x.*y;
         2*x.^2-x;
         2*y.^2-y;
         4*x-4*x.^2-4*x.*y;
         4*x.*y;
         4*y-4*y.^2-4*x.*y
        ];

elseif orderX==1&&orderY==0

    f = [-3+4*x+4*y;
         4*x-1;
         0;
         4-8*x-4*y;
         4*y;
         -4*y
        ];         

elseif orderX==0&&orderY==1

    f = [-3+4*y+4*x;
         0;
         4*y-1;
         -4*x;
         4*x;
         4-8*y-4*x
        ];

elseif orderX==2&&orderY==0  

    f = [4;
         4;
         0;
         -8;
         0;
         0
        ];

elseif orderX==0&&orderY==2 

    f = [4;
         0;
         4;
         0;
         0;
         -8
        ];

elseif orderX==1&&orderY==1 

    f = [4;
         0;
         0;
         -4;
         4;
         -4
        ];
      
end

end

function f = CrouzeixRaviart2D(x,y,orderX,orderY)

if orderX == 0 && orderY == 0
    
    f = [1-2*y; 2*x+2*y-1; 1-2*x];

elseif orderX == 1 && orderY == 0
    
    f = [0; 2; -2];

elseif orderX == 0 && orderY == 1
    
    f = [-2; 2; 0];
        
end

end