function A = skeletonAdjacencyMatrix(Skel)
% A = skeletonAdjacencyMatrix(Skel)
% Produces adjacency matrix of a skeleton.
% Variables:
% A - adjacency matrix.
% Skel - skeleton structure.
%
% David Pickup 2014

nJnts = numel(Skel.X);
I = Skel.E(:,1);
J = Skel.E(:,2);
S = ones(size(I));
A = sparse(I,J,S,nJnts,nJnts);
A = A + A';
A = double(A>0);

return;