function result = local_basis_2D(x,y,vertices,basis_type,basis_index,basis_der_x,basis_der_y)

% LOCAL_BASIS_2D
% 

x1 = vertices(1,1);
y1 = vertices(2,1);
x2 = vertices(1,2);
y2 = vertices(2,2);
x3 = vertices(1,3);
y3 = vertices(2,3);

J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

x_hat = ((y3-y1)*(x-x1)-(x3-x1)*(y-y1))/J;
y_hat = (-(y2-y1)*(x-x1)+(x2-x1)*(y-y1))/J;

if basis_type == 201
    
    psi_x = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,1,0);
    psi_y = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,1);
    
    if basis_der_x == 0 && basis_der_y == 0
        result = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,0);
    elseif basis_der_x == 1 && basis_der_y == 0
        result = psi_x*(y3-y1)/J+psi_y*(y1-y2)/J;
    elseif basis_der_x == 0 && basis_der_y == 1
        result = psi_x*(x1-x3)/J+psi_y*(x2-x1)/J;
    else
        result = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,basis_der_x,basis_der_y);
    end
    
elseif basis_type == 202
    
    psi_x = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,1,0);
    psi_y = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,1);
    psi_xx= reference_basis_2D(x_hat,y_hat,basis_type,basis_index,2,0);
    psi_xy= reference_basis_2D(x_hat,y_hat,basis_type,basis_index,1,1);
    psi_yy= reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,2);
    
    if basis_der_x == 0 && basis_der_y == 0
        result = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,0);
    elseif basis_der_x == 1 && basis_der_y == 0
        result = psi_x*(y3-y1)/J+psi_y*(y1-y2)/J;
    elseif basis_der_x == 0 && basis_der_y == 1
        result = psi_x*(x1-x3)/J+psi_y*(x2-x1)/J;
    elseif basis_der_x == 1 && basis_der_y == 1
        result = psi_xx*(x1-x3)*(y3-y1)/(J*J)+psi_xy*(x1-x3)*(y1-y2)/(J*J)+psi_xy*(x2-x1)*(y3-y1)/(J*J)+psi_yy*(x2-x1)*(y1-y2)/(J*J);
    elseif basis_der_x == 2 && basis_der_y == 0
        result = psi_xx*(y3-y1)*(y3-y1)/(J*J)+2*psi_xy*(y3-y1)*(y1-y2)/(J*J)+psi_yy*(y1-y2)*(y1-y2)/(J*J);
    elseif basis_der_x == 0 && basis_der_y == 2
        result = psi_xx*(x1-x3)*(x1-x3)/(J*J)+2*psi_xy*(x1-x3)*(x2-x1)/(J*J)+psi_yy*(x2-x1)*(x2-x1)/(J*J);
    else
        result = reference_basis_2D(x_hat,y_hat,basis_type,basis_index,basis_der_x,basis_der_y);
    end
    
end

end

