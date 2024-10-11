function int_value = Gauss_int_L_infi_error_1D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,basis_der_x)

Gauss_point_number = size(Gauss_weights,2);

int_value = 0;

for k = 1:Gauss_point_number
    
    value = analytic_solution(Gauss_nodes(k),basis_der_x)-local_FE_function_1D(number_of_local_basis_fun,uh_local_vec,Gauss_nodes(k),vertices,basis_type,basis_der_x);
    
    value = abs(value);
    
    if value > int_value
        int_value = value;
    end
    
end



end