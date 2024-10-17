function result = FE_basis_local_fun_1D(x,vertices,basis_type,basis_index,basis_der_x)
% 101:1D linear
% 102:1D quadratic

h = vertices(2) - vertices(1);

if basis_type == 101
    
    if basis_index == 1
        
        if basis_der_x == 0
            result = (vertices(2)-x)/h;
        elseif basis_der_x == 1
            result = -1/h;
        elseif basis_der_x >= 2 && (basis_der_x == fix(basis_der_x))
            result = 0;
        else
            warning('wrong input for basis derivative order!');
        end
        
    elseif basis_index == 2
        
        if basis_der_x == 0
            result = (x-vertices(1))/h;
        elseif basis_der_x == 1
            result = 1/h;
        elseif basis_der_x >= 2 && (basis_der_x == fix(basis_der_x))
            result = 0;
        else
            warning('wrong input for basis derivative order!');
        end
        
    else
        warning('wrong input for basis index!');
    end
    
% elseif basis_type == 102
%     
%     ?
    
end



end