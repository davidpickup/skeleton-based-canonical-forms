function Skel = AuSkeleton(mesh)
% Skel = AuSkeleton(mesh)
% Extracts a skeleton using the method by Au et al. SIGGRAPH 2008.
% Variables:
% Skel - extracted skeleton.
% mesh - input mesh structure.
%
% David Pickup 2014

% Perform mesh contraction.
out = AuMeshContraction(mesh);

% Produce 1D Skeleton.
A = adjacencyMatrix(out);
[I,J] = find(A);
[Skel.E, Skel.H] = AuConnectivitySurgeryMex(out.X,out.Y,out.Z,out.TRIV,J,I);
[idx, ~] = find(Skel.E ~= 0);
Skel.E = Skel.E(idx,:);
Skel.X = out.X;
Skel.Y = out.Y;
Skel.Z = out.Z;

% Delete duplicate edges.
Skel.E = sort(Skel.E,2);
Skel.E = unique(Skel.E,'rows');

% Delete unused vertices.
verts = unique(Skel.E);
newVerts = 1:numel(verts);
Skel.X = Skel.X(verts);
Skel.Y = Skel.Y(verts);
Skel.Z = Skel.Z(verts);
Skel.H = {Skel.H{verts}}';
for i = 1:numel(verts)
    Skel.E(Skel.E==verts(i)) = newVerts(i);
    Skel.H{i} = [verts(i); Skel.H{i}];
end

% Refine the skeleton's embedding.
Skel = AuEmbeddingRefinement(Skel, mesh);

return;

