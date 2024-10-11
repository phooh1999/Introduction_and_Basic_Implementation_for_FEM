function b = assemble_vector_2D(load_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_test,number_of_local_basis_fun_test,basis_type_test,basis_der_x_test,basis_der_y_test)

b = zeros(matrix_size(2),1);

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
%   T(:,n)第n个单元的所有点的索引
%   P(:,T(:,n))取出点的索引后，通过P矩阵取出具体坐标
%   在一维问题中也就是个1*2矩阵
    
    [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type);

    for beta = 1:number_of_local_basis_fun_test
            
            int_value = Gauss_quad_2D_test(load_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,beta,basis_der_x_test,basis_der_y_test);

            b(Tb_test(beta,n),1) = b(Tb_test(beta,n),1)+int_value;
%           beta是test function, alpha是trial function, 组装时需要反过来
    end


end



end