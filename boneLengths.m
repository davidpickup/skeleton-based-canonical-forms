function D = boneLengths(Skel)
% D = boneLengths(Skel)
% Produces adjacency matrix containing bone lengths.
% Variables:
% D - cost/adjacency matrix.
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

verts = [Skel.X Skel.Y Skel.Z];

[I,J] = find(A);
S = sqrt(sum((verts(I,:) - verts(J,:)).^2,2));
D = sparse(I,J,S,nJnts,nJnts);

return;