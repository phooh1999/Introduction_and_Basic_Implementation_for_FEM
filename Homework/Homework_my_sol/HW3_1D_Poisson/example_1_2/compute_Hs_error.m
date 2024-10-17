function [error] = compute_Hs_error(solution_vec,analytic_solution,P,T,Tb_trial,Gauss_type,number_of_local_basis_fun,basis_type,der_s)

% COMPUTE_HS_ERROR ¼ÆËãs½×H·¶
% s = 0,1

error = 0;
number_of_elements = size(T,2);

for n = 1:number_of_elements
    
    uh_local_vec = solution_vec(Tb_trial(:,n));
    
    vertices = P(:,T(:,n));
    
    [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type);
    
    int_value = Gauss_int_error_1D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,der_s);
    
    error = error + int_value;
    
end

error = sqrt(error);

end

