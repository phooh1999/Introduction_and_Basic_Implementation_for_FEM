function int_value = Gauss_int_error_2D(analytic_solution,Gauss_weights,Gauss_nodes,number_of_local_basis_fun,uh_local_vec,vertices,basis_type,basis_der_x,basis_der_y)

Gauss_point_number = size(Gauss_weights,2);

int_value = 0;

for k = 1:Gauss_point_number
    
    int_value = int_value + Gauss_weights(k)*(analytic_solution(Gauss_nodes(1,k),Gauss_nodes(2,k),basis_der_x,basis_der_y)-local_FE_function_2D(number_of_local_basis_fun,uh_local_vec,Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type,basis_der_x,basis_der_y))^2;
    
end



end