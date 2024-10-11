function [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type)
% 高斯求积公式局部
% local: 求出具体区间的求积点和权重
% Gauss_type == 10i: 1表示1维，i表示求积点个数
% Gauss_type == 20i: 2表示2维，i表示求积点个数
%% 一维
if Gauss_type == 102
    Gauss_weights = [1,1];
    Gauss_nodes = [-1/sqrt(3),1/sqrt(3)];
elseif Gauss_type == 103
    Gauss_weights = [5/9,8/9,5/9];
    Gauss_nodes = [-sqrt(0.6),0,sqrt(0.6)];
elseif Gauss_type == 104
    Gauss_weights = [0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454];
    Gauss_nodes = [-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053];
%% 二维
elseif Gauss_type == 209

    Gauss_ref_weights=[64/81*(1-0)/8,100/324*(1-sqrt(3/5))/8,100/324*(1-sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,40/81*(1-0)/8,40/81*(1-0)/8,40/81*(1-sqrt(3/5))/8,40/81*(1+sqrt(3/5))/8];
    Gauss_ref_nodes=[(1+0)/2,(1-0)*(1+0)/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+sqrt(3/5))/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1-sqrt(3/5))/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+sqrt(3/5))/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1-sqrt(3/5))/4;(1+0)/2,(1-0)*(1+sqrt(3/5))/4;(1+0)/2,(1-0)*(1-sqrt(3/5))/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+0)/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+0)/4];

%%
else
    warning('Wrong Gauss Type!');
end

if Gauss_type == 209
    
    x1=vertices(1,1);
    y1=vertices(2,1);
    x2=vertices(1,2);
    y2=vertices(2,2);
    x3=vertices(1,3);
    y3=vertices(2,3);
    
    Jacobi=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    
    Gauss_weights = Jacobi*Gauss_ref_weights;
    Gauss_nodes(:,1) = x1+(x2-x1)*Gauss_ref_nodes(:,1)+(x3-x1)*Gauss_ref_nodes(:,2);
    Gauss_nodes(:,2) = y1+(y2-y1)*Gauss_ref_nodes(:,1)+(y3-y1)*Gauss_ref_nodes(:,2);
    Gauss_nodes = Gauss_nodes';

else
    Gauss_weights(1,:) = (vertices(2) - vertices(1))/2 * Gauss_weights(1,:);
    Gauss_nodes(1,:) = (vertices(2)+vertices(1))/2 + (vertices(2)-vertices(1))/2*Gauss_nodes(:);

end

end