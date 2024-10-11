clear
clc

%%

x = [0,1,0,1];
N = [4,4];

mesh_type = 202;

[P,T] = generate_PT(x,N,201);
[Pb,Tb] = generate_PT(x,N,mesh_type);
[boundaryedges] = generate_boundaryedges(N,T);
[boundaryedges1] = generate_boundaryedges(N,Tb);
[boundarynodes] = generate_boundarynodes(boundaryedges1,mesh_type);
vertices = P(:,T(:,3));

%% check

[boundary_nodes_1,boundary_edges_1] = generate_boundary_nodes_edges(2*N(1),2*N(2),N(1),N(2));
boundary_nodes = boundarynodes - boundary_nodes_1;
boundary_edges = boundaryedges - boundary_edges_1;

[Mc,Tc]=generate_M_T_triangle(0,1,0,1,[0.25,0.25],2);

Pc = Pb - Mc;
Tcc = Tb - Tc;

