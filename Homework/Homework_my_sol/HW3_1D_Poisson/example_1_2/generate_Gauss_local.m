function [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type)
% 高斯求积公式局部
% local: 求出具体区间的求积点和权重
% Gauss_type == 10i: 1表示1维，i表示求积点个数

if Gauss_type == 102
    Gauss_weights = [1,1];
    Gauss_nodes = [-1/sqrt(3),1/sqrt(3)];
elseif Gauss_type == 103
    Gauss_weights = [5/9,8/9,5/9];
    Gauss_nodes = [-sqrt(0.6),0,sqrt(0.6)];
elseif Gauss_type == 104
    Gauss_weights = [0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454];
    Gauss_nodes = [-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053];
%     Gauss_weights = [0.3478548451,0.6521451549,0.6521451549,0.3478548451];
%     Gauss_nodes = [-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
else
    warning('Wrong Gauss Type!');
end

Gauss_weights(1,:) = (vertices(2) - vertices(1))/2 * Gauss_weights;
Gauss_nodes(1,:) = (vertices(2)+vertices(1))/2 + (vertices(2)-vertices(1))/2*Gauss_nodes(:);

end