function mesh = skeletonCanonicalForm(mesh)
% mesh = skeletonCanonicalForm(mesh)
% Computes a canonical form of a mesh using its skeleton.
%
% David Pickup 2015

mesh = normaliseMesh(mesh,100);

% MDS options.
options.rtol = 0.001;
options.iter = 100;
options.method = 'smacof';
options.dim = 3;

% Extract skeleton from mesh.
Skel = AuSkeleton(mesh);
Skel.TBA = triangleBoneAssignment(Skel, mesh);

% Make sure the skeleton is fully connected.
 D = boneLengths(Skel);
D = graphallshortestpaths(D);
if numel(find(isinf(D))) > 0
    while numel(find(isinf(D))) > 0
        DE = allPairsEuclideanMesh(Skel);
        DE(~isinf(D)) = Inf;
        [I,J] = find(DE==min(DE(:)));
        Skel.E(end+1,:) = [I(1) J(1)];
        D = boneLengths(Skel);
        D = allPairsDijkstra(D);
    end
    %Skel = mergeJunctions(Skel,mesh);
    Skel = AuEmbeddingRefinement(Skel, mesh);
    Skel.TBA = triangleBoneAssignment(Skel, mesh);
    D = boneLengths(Skel);
    D = graphallshortestpaths(D);
end

% Perform MDS on the skeleton.
verts = [Skel.X Skel.Y Skel.Z];
options.X0 = verts;
X = mds(D,options);

% Compute the mesh deformation based on the transformed skeleton.
Skel2 = Skel;
Skel2.X = X(:,1);
Skel2.Y = X(:,2);
Skel2.Z = X(:,3);
mesh = YanSkeletonDrivenMeshDeformation(mesh, Skel, Skel2);

return;