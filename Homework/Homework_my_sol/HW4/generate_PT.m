function [P,T] = generate_PT(x,N,mesh_type)
% x: ��������
% N: ����������
% mesh_type: ���񻮷�����
% mesh_type == 101: 1D uniform partition
%% һά����
if mesh_type == 101
    P = zeros(1,N+1);
    h = (x(2) - x(1))/N;
    P(1,:) = x(1):h:x(2);
    
    T = zeros(2,N);
    T(1,:) = 1:N;
    T(2,:) = 2:(N+1);
%% ��Ԫ���е�
elseif mesh_type == 102
    P = zeros(1,2*N+1);
    h = (x(2) - x(1))/(2*N);
    P(1,:) = x(1):h:x(2);
    
    T = zeros(3,N);
    T(1,:) = 1:2:(2*N-1);
    T(2,:) = 3:2:(2*N+1);
    T(3,:) = 2:2:2*N;
%   Tb(1,k)��k����Ԫ����˵�
%   Tb(2,k)��k����Ԫ���Ҷ˵�
%   Tb(3,k)��k����Ԫ���е�
%% ��ά����������
elseif mesh_type == 201
    
    r_n = N(2)+1;
    c_n = N(1)+1;
    
    number_of_nodes = r_n*c_n;
    number_of_elements = 2*N(1)*N(2);
    
    P = zeros(2,number_of_nodes);
    T = zeros(3,number_of_elements);
    h = [(x(2)-x(1))/N(1),(x(4)-x(3))/N(2)];
    
    for i = 1:r_n
        for j = 1:c_n
            x0 = x(1)+(j-1)*h(1);
            y = x(3)+(i-1)*h(2);
            n = (j-1)*r_n+i;
            P(1,n) = x0;
            P(2,n) = y;
        end
    end
    
    for i = 1:N(1)
        for j = 1:N(2)
            index = (i-1)*2*N(2) + (j-1)*2;
            T(:,index+1) = [(i-1)*r_n+j;i*r_n+j;(i-1)*r_n+j+1];
            T(:,index+2) = [(i-1)*r_n+j+1;i*r_n+j;i*r_n+j+1];
        end
    end
    
%% ��ά������������е㣬��������Ԫ����
elseif mesh_type == 202
    
    r_n = 2*N(2)+1;
    c_n = 2*N(1)+1;
    
    number_of_nodes = r_n*c_n;
    number_of_elements = 2*N(1)*N(2);
    
    P = zeros(2,number_of_nodes);
    T = zeros(6,number_of_elements);
    h = [(x(2)-x(1))/(2*N(1)),(x(4)-x(3))/(2*N(2))];
    
    for i = 1:r_n
        for j = 1:c_n
            x0 = x(1)+(j-1)*h(1);
            y = x(3)+(i-1)*h(2);
            n = (j-1)*r_n+i;
            P(1,n) = x0;
            P(2,n) = y;
        end
    end
    
    for i = 1:N(1)
        for j = 1:N(2)
            index = (i-1)*2*N(2) + (j-1)*2;
            vertice = 2*(i-1)*r_n+2*(j-1)+1;
            T(:,index+1) = [vertice;vertice+2*r_n;vertice+2;vertice+r_n;vertice+r_n+1;vertice+1];
            T(:,index+2) = [vertice+2;vertice+2*r_n;vertice+2*r_n+2;vertice+r_n+1;vertice+2*r_n+1;vertice+r_n+2];
        end
    end
    
end

end

