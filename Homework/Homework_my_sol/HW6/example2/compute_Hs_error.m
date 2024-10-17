function [error] = compute_Hs_error(solution_vec,analytic_solution,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun,basis_type,der_s)

% COMPUTE_HS_ERROR ¼ÆËãs½×H·¶
% s = 0,1

error = 0;
number_of_elements = size(T,2);

for n = 1:number_of_elements
    
    uh_local_vec = solution_vec(Tb_trial(:,n));
    
    vertices = P(:,T(:,n));
    
    [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type);
    
    if der_s == 0
    
        int_value = Gauss_int_error_2D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,der_s,der_s);
    
        error = error + int_value;
    
    elseif der_s == 1
        
        int_value1 = Gauss_int_error_2D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,der_s,0);
        int_value2 = Gauss_int_error_2D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,0,der_s);
        
        error = error + int_value1 + int_value2;
    
    end
    
end

error = sqrt(error);

end

