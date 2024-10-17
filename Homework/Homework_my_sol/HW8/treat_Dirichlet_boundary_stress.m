function [A,b] = treat_Dirichlet_boundary_stress(boundarynodes,bound_fun,Pb,A,b)

nbn = size(boundarynodes,2);
nb = size(Pb,2);

for k = 1:nbn
    
    if boundarynodes(1,k) == -1
        
        i = boundarynodes(2,k);
        A(i,:) = 0;
        A(i,i) = 1;
        x = Pb(1,i);
        y = Pb(2,i);
        A(nb+i,:) = 0;
        A(nb+i,nb+i) = 1;
        [b(i),b(nb+i)] = bound_fun(x,y);
        
    end
    
end


end