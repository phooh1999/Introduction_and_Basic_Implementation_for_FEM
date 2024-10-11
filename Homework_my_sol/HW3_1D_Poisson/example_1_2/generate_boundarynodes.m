function boundarynodes = generate_boundarynodes(x_domain,mesh_type)

if mesh_type == 101
    boundarynodes(:,1) = x_domain(1);
    boundarynodes(:,2) = x_domain(2);
end

end