function TriMeshGenerateTest(type)

problemDomain.x = [0,1];
problemDomain.y = [0,1];

partitionSteps.x = 10;
partitionSteps.y = 8;

mesh = TriMeshGenerate(problemDomain, partitionSteps, type);

hold on

trimesh(mesh.T(:,1:3), mesh.P(:,1), mesh.P(:,2));
scatter(mesh.P(:,1), mesh.P(:,2));

x = mesh.P(:,1);
y = mesh.P(:,2);

for i = 1:max(size(x))
    c = num2str(i-1);
    text(x(i), y(i),c);
end

end