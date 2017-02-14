function mesh = normaliseMesh(mesh, s)
% mesh = normaliseMesh(mesh, [scale])
% Centre mesh and scale to unit size.
%
% David Pickup 2013

% Translate mesh so its centroid is at the origin.
mesh.X = mesh.X - mean(mesh.X);
mesh.Y = mesh.Y - mean(mesh.Y);
mesh.Z = mesh.Z - mean(mesh.Z);

% Scale mesh so the largest radial distance is equal to 1.
distance = sqrt(mesh.X.^2 + mesh.Y.^2 + mesh.Z.^2);
scale = 1/max(distance);
mesh.X = mesh.X .* scale;
mesh.Y = mesh.Y .* scale;
mesh.Z = mesh.Z .* scale;
clear distance;

if nargin == 2
    mesh.X = mesh.X .* s;
    mesh.Y = mesh.Y .* s;
    mesh.Z = mesh.Z .* s;
end