function result = FE_basis_local_fun_1D(x,vertices,basis_type,basis_index,basis_der_x)
% 101:1D linear
% 102:1D quadratic

h = vertices(2) - vertices(1);

%% 一阶单元
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
    
%% 二阶单元
elseif basis_type == 102
    
    if basis_index == 1
        
        if basis_der_x == 0
            result = 2*((x-vertices(1))/h)^2-3*((x-vertices(1))/h)+1;
        elseif basis_der_x == 1
            result = 4*((x-vertices(1))/h)*1/h-3/h;
        elseif basis_der_x == 2 
            result = 4/(h*h);
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        else
            warning('wrong input for basis derivative order!');
        end
        
    elseif basis_index == 2
        
        if basis_der_x == 0
            result = 2*((x-vertices(1))/h)^2-((x-vertices(1))/h);
        elseif basis_der_x == 1
            result = 4*((x-vertices(1))/h)*1/h-1/h;
        elseif basis_der_x == 2
            result = 4/(h*h);
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        else
            warning('wrong input for basis derivative order!');
        end
        
    elseif basis_index == 3
        
        if basis_der_x == 0
            result = -4*((x-vertices(1))/h)^2+4*((x-vertices(1))/h);
        elseif basis_der_x == 1
            result = -8*((x-vertices(1))/h)*1/h+4/h;
        elseif basis_der_x == 2
            result = -8/(h*h);
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        else
            warning('wrong input for basis derivative order!');
        end
        
    else
        warning('wrong input for basis index!');
    end
    
end



end