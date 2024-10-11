function [A,b] = treat_Dirichlet_boundary(boundarynodes,bound_fun,Pb,A,b)

nbn = size(boundarynodes,2);

for k = 1:nbn
    
    if boundarynodes(1,k) == -1
        
        i = boundarynodes(2,k);
        A(i,:) = 0;
        A(i,i) = 1;
        x = Pb(1,i);
        y = Pb(2,i);
        b(i) = bound_fun(x,y);
        
    end
    
end


end