function result = reference_basis_2D(x,y,basis_type,basis_index,basis_der_x,basis_der_y)

% REFERENCE_BASIS_2D 
% 二维单元参考基函数
% 201: 2D linear
% 202: 2D quadratic

%% 2D linear
if basis_type == 201
    
    if basis_index == 1
        
        if basis_der_x == 0 && basis_der_y == 0
            result = -x-y+1;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = -1;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = -1;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 0;
        elseif basis_der_x >= 2 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 2 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 2
        
        if basis_der_x == 0 && basis_der_y == 0
            result = x;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 1;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 0;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 0;
        elseif basis_der_x >= 2 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 2 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 3
        
        if basis_der_x == 0 && basis_der_y == 0
            result = y;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 0;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 1;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 0;
        elseif basis_der_x >= 2 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 2 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    else
        warning('Wrong Input for basis_index!')
    end
    
%% 2D quadratic
elseif basis_type == 202
    
    if basis_index == 1
        
        if basis_der_x == 0 && basis_der_y == 0
            result = 2*x*x+2*y*y+4*x*y-3*x-3*y+1;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 4*x+4*y-3;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 4*x+4*y-3;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 4;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = 4;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = 4;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 2
        
        if basis_der_x == 0 && basis_der_y == 0
            result = 2*x*x-x;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 4*x-1;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 0;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 0;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = 4;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = 0;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 3
        
        if basis_der_x == 0 && basis_der_y == 0
            result = 2*y*y-y;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 0;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 4*y-1;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 0;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = 4;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 4
        
        if basis_der_x == 0 && basis_der_y == 0
            result = -4*x*x-4*x*y+4*x;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = -8*x-4*y+4;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = -4*x;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = -4;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = -8;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = 0;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
            
    elseif basis_index == 5
        
        if basis_der_x == 0 && basis_der_y == 0
            result = 4*x*y;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = 4*y;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = 4*x;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = 4;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = 0;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    elseif basis_index == 6
        
        if basis_der_x == 0 && basis_der_y == 0
            result = -4*y*y-4*x*y+4*y;
        elseif basis_der_x == 1 && basis_der_y == 0
            result = -4*y;
        elseif basis_der_x == 0 && basis_der_y == 1
            result = -8*y-4*x+4;
        elseif basis_der_x == 1 && basis_der_y == 1
            result = -4;
        elseif basis_der_x == 2 && basis_der_y == 0
            result = 0;
        elseif basis_der_x == 0 && basis_der_y == 2
            result = -8;
        elseif basis_der_x >= 3 && (basis_der_x == fix(basis_der_x))
            result = 0;
        elseif basis_der_y >= 3 && (basis_der_y == fix(basis_der_y))
            result = 0;
        else
            warning('Wrong Input for basis_der!')
        end
        
    else
        
        warning('Wrong Input for basis_index!')
            
    end    
    
end

end

