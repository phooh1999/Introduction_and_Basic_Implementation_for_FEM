function v = treat_Neumann_boundary(size_of_v,boundaryedges,load_fun,P,T,Tb_test,number_of_local_basis_fun,basis_type_test)

nbn = size(boundaryedges,2);

v = zeros(size_of_v,1);



for k = 1:nbn
    
    if boundaryedges(1,k) == -2
        
        n_k = boundaryedges(2,k);
        
        edge = P(:,boundaryedges(3:4,k));
        
        vertices = P(:,T(:,n_k));
        
        [Gauss_weights,Gauss_nodes] = generate_Gauss_local_1D_line(edge,104);
        
        for beta = 1:number_of_local_basis_fun
            
            int_value = Gauss_quad_2D_test(load_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,beta,0,0);
            
            v(Tb_test(beta,n_k),1) = v(Tb_test(beta,n_k),1) + int_value;
            
        end
        
    end
    
end



end