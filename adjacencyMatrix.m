function A = adjacencyMatrix(mesh)
% A = adjacencyMatrix(mesh)
% Computes the adjacency matrix of a mesh.
% Variables:
% A - adjacency matrix.
% mesh - mesh structure.
%
% David Pickup 2013

% Get number of vertices.
nPnts = numel(mesh.X);

% Compute adjacency matrix.
I = [mesh.TRIV(:,1);mesh.TRIV(:,2);mesh.TRIV(:,3)];
J = [mesh.TRIV(:,2);mesh.TRIV(:,3);mesh.TRIV(:,1)];
S = ones(size(I));
A = sparse(I,J,S,nPnts,nPnts);
A = A + A';
A = double(A>0);

return;