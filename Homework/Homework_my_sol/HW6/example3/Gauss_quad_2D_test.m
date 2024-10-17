function int_value = Gauss_quad_2D_test(coe_fun,Gauss_weights,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_test,basis_der_y_test)

Gauss_point_number = size(Gauss_weights,2);

int_value = 0;

for k = 1:Gauss_point_number
    
    int_value = int_value + Gauss_weights(1,k)*feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_test,basis_index_test,basis_der_x_test,basis_der_y_test);
    
end



end