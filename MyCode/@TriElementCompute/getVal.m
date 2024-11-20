function val = getVal(obj,x,y,orderX,orderY,femSol,vertices)

basis = obj.trialFun.evaluate(x,y,vertices,orderX,orderY);

val = femSol' * basis;

end

