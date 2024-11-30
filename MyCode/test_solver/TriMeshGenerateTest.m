function TriMeshGenerateTest(type)

problemDomain.x = [0,1];
problemDomain.y = [0,1];

partitionSteps.x = 4;
partitionSteps.y = 5;

mesh = TriMeshGenerate(problemDomain, partitionSteps, 1);
element = TriMeshGenerate(problemDomain, partitionSteps, type);

boundary = TriMeshBoundary(partitionSteps, type);

hold on
axis off

x = element.P(:,1);
y = element.P(:,2);

trimesh(mesh.T(:,1:3), mesh.P(:,1), mesh.P(:,2));
scatter(x,y);

for i = 1:max(size(x))
    c = num2str(i);
    text(x(i), y(i),c);
end

end