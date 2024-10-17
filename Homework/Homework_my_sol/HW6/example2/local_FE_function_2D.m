function result = local_FE_function_2D(number_of_local_basis_fun,uh_local_vec,x,y,vertices,basis_type,basis_der_x,basis_der_y)

result = 0;

for k = 1:number_of_local_basis_fun
    
    result = result + uh_local_vec(k)*local_basis_2D(x,y,vertices,basis_type,k,basis_der_x,basis_der_y);
    
end


end