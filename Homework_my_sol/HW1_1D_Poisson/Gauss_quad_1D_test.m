function int_value = Gauss_quad_1D_test(coe_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_test)

Gauss_point_number = size(Gauss_weights,2);

int_value = 0;

for k = 1:Gauss_point_number
    
    int_value = int_value + Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k))*FE_basis_local_fun_1D(Gauss_nodes(k),vertices,basis_type_test,basis_index_test,basis_der_x_test);
    
end



end