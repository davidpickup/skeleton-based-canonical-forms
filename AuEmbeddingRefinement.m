function Skel = AuEmbeddingRefinement(Skel, mesh)
% Skel = AuEmbeddingRefinement(Skel, mesh)
% Refines the embedding of a skeleton using the method by
% Au et al. SIGGRAPH 2008.
% Variables:
% Skel - skeleton structure.
% mesh - mesh structure.
%
% David Pickup 2014

% Skel = mergeJunctions(Skel,mesh);

nPnts = numel(Skel.X);

Adj = adjacencyMatrix(mesh);
[MJ,MI] = find(Adj);

%% Merge junctions.
flag = 1;
while flag
    oldE = Skel.E;
    skelAdj = skeletonAdjacencyMatrix(Skel);

    [SJ,SI] = find(skelAdj);
    [Skel.E,Skel.H] = AuMergeJunctionsMex(mesh.X,mesh.Y,mesh.Z,Skel.X,Skel.Y,Skel.Z,SI,SJ,Skel.H);

    % Delete duplicate edges.
    [idx, ~] = find(Skel.E ~= 0);
    Skel.E = Skel.E(idx,:);
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
    end
    
    if numel(Skel.E) == numel(oldE)
        flag = 0;
    end
end

%% Centralise skeleton.
skelAdj = skeletonAdjacencyMatrix(Skel);
[SJ,SI] = find(skelAdj);
[X,Y,Z] = AuEmbeddingRefinementMex(mesh.X,mesh.Y,mesh.Z,Skel.X,Skel.Y,Skel.Z,MI,MJ,SI,SJ,Skel.H);
idx = find(isnan(X+Y+Z));
X(idx) = Skel.X(idx);
Y(idx) = Skel.Y(idx);
Z(idx) = Skel.Z(idx);
Skel.X = X;
Skel.Y = Y;
Skel.Z = Z;







return;