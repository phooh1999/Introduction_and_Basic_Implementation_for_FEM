function A = assemble_matrix_1D(coe_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test)

A = sparse(matrix_size(1),matrix_size(2));

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
%   T(:,n)第n个单元的所有点的索引
%   P(:,T(:,n))取出点的索引后，通过P矩阵取出具体坐标
%   在一维问题中也就是个1*2矩阵
    
    [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type);
    
    for alpha = 1:number_of_local_basis_fun_trial
        
        for beta = 1:number_of_local_basis_fun_test
            
            int_value = Gauss_quad_1D_trial_test(coe_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_trial,alpha,basis_der_x_trial,basis_type_test,beta,basis_der_x_test);
            
%             S(beta,aloha) = int_value;
%             之后可以改成单刚组装总刚

            A(Tb_test(beta,n),Tb_trial(alpha,n)) = A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
%           beta是test function, alpha是trial function, 组装时需要反过来
        end
        
    end
    
end




end