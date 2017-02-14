function L = laplaceBeltrami(mesh)
% L = laplaceBeltrami(mesh)
% Computes the Laplace-Beltrami matrix using cotangent weights. Uses code
% by Ben-Chen.
% Variables:
% L - Laplace-Beltrami matrix.
% mesh - mesh structure.
%
% David Pickup 2013

T = mesh.TRIV;
nv = numel(mesh.X);
X = [mesh.X,mesh.Y,mesh.Z];

L1 = len(X(T(:,2),:)-X(T(:,3),:));
L2 = len(X(T(:,1),:)-X(T(:,3),:));
L3 = len(X(T(:,1),:)-X(T(:,2),:));
EL = [L1, L2, L3];

epsi = 1e-8;
L1(find(L1==0)) = epsi;
L2(find(L2==0)) = epsi;
L3(find(L3==0)) = epsi;

A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);

A = [A1,A2,A3];
A = acos(A);

A = real(A);
A(find(A==0)) = epsi;

I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);

W = sparse(I,J,S,nv,nv);

nv = length(W);

[I,J,S] = find(W);
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [-S;-S;S;S];

L = sparse(In,Jn,Sn,nv,nv);

return;

function n = len(X)
n = sqrt(sum(X.^2,2));
