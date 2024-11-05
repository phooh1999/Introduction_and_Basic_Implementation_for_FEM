function res = evaluate(obj,x,y,vertices,orderX,orderY)

J_11=vertices(2,1)-vertices(1,1);
J_12=vertices(3,1)-vertices(1,1);
J_21=vertices(2,2)-vertices(1,2);
J_22=vertices(3,2)-vertices(1,2);
J_det=J_11*J_22-J_12*J_21;

x_hat=(J_22*(x-vertices(1,1))-J_12*(y-vertices(1,2)))/J_det;
y_hat=(-J_21*(x-vertices(1,1))+J_11*(y-vertices(1,2)))/J_det;

if orderX==0&&orderY==0
    res=obj.refBasisFun(x_hat,y_hat,0,0);
elseif orderX==1&&orderY==0
    res=(obj.refBasisFun(x_hat,y_hat,1,0).*J_22 ...
       + obj.refBasisFun(x_hat,y_hat,0,1).*(-J_21))./J_det;
elseif orderX==0&&orderY==1
    res=(obj.refBasisFun(x_hat,y_hat,1,0).*(-J_12) ...
       + obj.refBasisFun(x_hat,y_hat,0,1).*J_11)./J_det;
elseif orderX==2&&orderY==0
    res=(obj.refBasisFun(x_hat,y_hat,2,0).*J_22^2 ...
       + obj.refBasisFun(x_hat,y_hat,0,2).*J_21^2 ...
       + obj.refBasisFun(x_hat,y_hat,1,1).*(-2*J_21*J_22))./J_det^2;
elseif orderX==0&&orderY==2
    res=(obj.refBasisFun(x_hat,y_hat,2,0).*J_12^2 ...
       + obj.refBasisFun(x_hat,y_hat,0,2).*J_11^2 ...
       + obj.refBasisFun(x_hat,y_hat,1,1).*(-2*J_11*J_12))./J_det^2;
elseif orderX==1&&orderY==1
    res=(obj.refBasisFun(x_hat,y_hat,2,0).*(-J_22*J_12) ...
       + obj.refBasisFun(x_hat,y_hat,0,2).*(-J_21*J_11) ...
       + obj.refBasisFun(x_hat,y_hat,1,1).*(J_21*J_12+J_11*J_22))./J_det^2;
end

end

