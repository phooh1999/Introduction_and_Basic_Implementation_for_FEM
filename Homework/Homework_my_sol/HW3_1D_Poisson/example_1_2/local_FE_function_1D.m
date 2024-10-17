function result = local_FE_function_1D(number_of_local_basis_fun,uh_local_vec,x,vertices,basis_type,basis_der_x)

result = 0;

for k = 1:number_of_local_basis_fun
    
    result = result + uh_local_vec(k)*FE_basis_local_fun_1D(x,vertices,basis_type,k,basis_der_x);
    
end


end