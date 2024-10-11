function u = generate_initial_value(Pb,initial_fun)

nbn = size(Pb,2);

u = zeros(nbn,1);

for i = 1:nbn
    
    x = Pb(1,i);
    y = Pb(2,i);
    
    u(i,1) = initial_fun(x,y);
    
end

end