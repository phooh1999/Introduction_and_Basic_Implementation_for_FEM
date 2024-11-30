function generatePT(obj, domain, steps, type)
% 1: 2D Lagrange linear FE
% 2: 2D Lagrange quadratic FE
% 3: 2D Crouzeix Raviart FE

%% 2D Lagrange linear FE
if type == 1
    x = [domain.x, domain.y];
    N = [steps.x, steps.y];
    
    r_n = N(2)+1;
    c_n = N(1)+1;
    
    number_of_nodes = r_n*c_n;
    number_of_elements = 2*N(1)*N(2);
    
    obj.P = zeros(2,number_of_nodes);
    obj.T = zeros(3,number_of_elements);
    h = [(x(2)-x(1))/N(1),(x(4)-x(3))/N(2)];
    
    for i = 1:r_n
        for j = 1:c_n
            x0 = x(1)+(j-1)*h(1);
            y = x(3)+(i-1)*h(2);
            n = (j-1)*r_n+i;
            obj.P(1,n) = x0;
            obj.P(2,n) = y;
        end
    end
    
    for i = 1:N(1)
        for j = 1:N(2)
            index = (i-1)*2*N(2) + (j-1)*2;
            obj.T(:,index+1) = [(i-1)*r_n+j;i*r_n+j;(i-1)*r_n+j+1];
            obj.T(:,index+2) = [(i-1)*r_n+j+1;i*r_n+j;i*r_n+j+1];
        end
    end
    
    obj.P = obj.P';
    obj.T = obj.T';
    return
end

%% 2D Lagrange quadratic FE
if type == 2
    x = [domain.x, domain.y];
    N = [steps.x, steps.y];
    
    r_n = 2*N(2)+1;
    c_n = 2*N(1)+1;
    
    number_of_nodes = r_n*c_n;
    number_of_elements = 2*N(1)*N(2);
    
    obj.P = zeros(2,number_of_nodes);
    obj.T = zeros(6,number_of_elements);
    h = [(x(2)-x(1))/(2*N(1)),(x(4)-x(3))/(2*N(2))];
    
    for i = 1:r_n
        for j = 1:c_n
            x0 = x(1)+(j-1)*h(1);
            y = x(3)+(i-1)*h(2);
            n = (j-1)*r_n+i;
            obj.P(1,n) = x0;
            obj.P(2,n) = y;
        end
    end
    
    for i = 1:N(1)
        for j = 1:N(2)
            index = (i-1)*2*N(2) + (j-1)*2;
            vertice = 2*(i-1)*r_n+2*(j-1)+1;
            obj.T(:,index+1) = [vertice;vertice+2*r_n;vertice+2;vertice+r_n;vertice+r_n+1;vertice+1];
            obj.T(:,index+2) = [vertice+2;vertice+2*r_n;vertice+2*r_n+2;vertice+r_n+1;vertice+2*r_n+1;vertice+r_n+2];
        end
    end
    
    obj.P = obj.P';
    obj.T = obj.T';
    return
end

%% 2D Crouzeix Raviart FE

if type == 3
    N1 = steps.x;        N2 = steps.y;
    x0 = domain.x(1);    x1 = domain.x(2);
    y0 = domain.y(1);    y1 = domain.y(2);
    dx = (x1-x0)/N1;     dy = (y1-y0)/N2;
    
    number_of_nodes = N2*(N1+1) + N1*(2*N2+1);
    number_of_elements = 2*N1*N2;
    
    obj.P = zeros(number_of_nodes,2);
    obj.T = zeros(number_of_elements,3);
    
    for i = 0 : N1
        for j = 1 : N2
            obj.P(i*(3*N2+1)+j,1) = x0 + i*dx;
            obj.P(i*(3*N2+1)+j,2) = y0 + (j-1/2)*dy;
        end
    end
    
    for i = 1 : N1
        for k = 1 : (2*N2+1)
            obj.P((3*i-2)*N2+i-1+k,1) = x0 + (i-1/2)*dx;
            obj.P((3*i-2)*N2+i-1+k,2) = y0 + (k-1)*1/2*dy;
        end
    end
    
    for i = 1 : N1
        for j = 1 : N2
            obj.T((i-1)*2*N2+2*j-1,1) = (i-1)*(3*N2+1) + N2 + 2*j - 1;
            obj.T((i-1)*2*N2+2*j-1,2) = (i-1)*(3*N2+1) + N2 + 2*j;
            obj.T((i-1)*2*N2+2*j-1,3) = (i-1)*(3*N2+1) + j;
            
            obj.T((i-1)*2*N2+2*j,1) = (i-1)*(3*N2+1) + N2 + 2*j;
            obj.T((i-1)*2*N2+2*j,2) = i*(3*N2+1) + j;
            obj.T((i-1)*2*N2+2*j,3) = (i-1)*(3*N2+1) + N2 + 2*j + 1;
        end
    end
    
    return
end

%% WARNING!!!

warning('Wrong TYPE of mesh!');

end

