function [nodes,jac] = localNodes(obj, vertices)

x1=vertices(1,1);   y1=vertices(1,2);
x2=vertices(2,1);   y2=vertices(2,2);
x3=vertices(3,1);   y3=vertices(3,2);
            
jac=0.5 * abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

nodes(:,1) = x1+(x2-x1)*obj.lroots(:,1)+(x3-x1)*obj.lroots(:,2);
nodes(:,2) = y1+(y2-y1)*obj.lroots(:,1)+(y3-y1)*obj.lroots(:,2);

end

