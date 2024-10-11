function [R,w] = treat_Robin_boundary(matrix_size,boundaryedges,robin_r,robin_q,P,T,Tb_test,Tb_trial,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test)

nbn = size(boundaryedges,2);

R = sparse(matrix_size(1,1),matrix_size(1,2));
w = zeros(matrix_size(1,2),1);

for k = 1:nbn
    
    if boundaryedges(1,k) == -3
        
        n_k = boundaryedges(2,k);
        
        edge = P(:,boundaryedges(3:4,k));
        
        vertices = P(:,T(:,n_k));
        
        [Gauss_weights,Gauss_nodes] = generate_Gauss_local_1D_line(edge,104);
        
        for beta = 1:number_of_local_basis_fun_test
            
            int_value = Gauss_quad_2D_test(robin_q,Gauss_weights,Gauss_nodes,vertices,basis_type_test,beta,0,0);
            
            w(Tb_test(beta,n_k),1) = w(Tb_test(beta,n_k),1) + int_value;
            
        end
        
        for alpha = 1:number_of_local_basis_fun_trial
            
            for beta = 1:number_of_local_basis_fun_test
                
                int_value = Gauss_quad_2D_trial_test(robin_r,Gauss_weights,Gauss_nodes,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
                
                R(Tb_test(beta,n_k),Tb_trial(alpha,n_k)) = R(Tb_test(beta,n_k),Tb_trial(alpha,n_k)) + int_value;
                
            end
            
        end
        
    end
    
end



end