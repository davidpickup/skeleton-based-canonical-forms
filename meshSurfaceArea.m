function [A,AF,AV] = meshSurfaceArea(mesh)
% [A,AF,AV] = meshSurfaceArea(mesh)
% Calculates the mesh surface area.
% Variables:
% A - total mesh surface area.
% AF - surface area of each face.
% AV - surface area of each vertex.
% mesh - mesh structure.
%
% David Pickup 2013

nPnts = numel(mesh.X);

verts = [mesh.X mesh.Y mesh.Z];

% Calculate the surface area per face.
E12 = verts(mesh.TRIV(:,1),:)-verts(mesh.TRIV(:,2),:);
E13 = verts(mesh.TRIV(:,1),:)-verts(mesh.TRIV(:,3),:);
N = cross(E12,E13);
AF = 0.5*normv(N);

% Calculate the surface area per vertex.
I = [mesh.TRIV(:,1);mesh.TRIV(:,2);mesh.TRIV(:,3)];
J = ones(length(I),1);
S = [AF;AF;AF];
AV = sparse(I,J,S,nPnts,1);

% Calculate the total surface area.
A = sum(AF);

return;

function nn = normv(V)
nn = sqrt(sum(V.^2,2));