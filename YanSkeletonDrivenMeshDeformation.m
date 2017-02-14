function mesh = YanSkeletonDrivenMeshDeformation(mesh, Skel1, Skel2)
% mesh = YanSkeletonDrivenMeshDeformation(mesh, Skel1, Skel2)
% Deforms a mesh based on a skeleton transformation using the method by Yan
% et al. TVCG 2008.
% Variables:
% mesh - mesh structure.
% Skel1 - initial skeleton structure.
% Skel2 - transformed skeleton structure.
%
% David Pickup 2014

nTris = size(mesh.TRIV,1);
nPnts = numel(mesh.X);
nOrigPnts = nPnts;
verts = [mesh.X mesh.Y mesh.Z]';

% Get surface area of each triangle.
[~,A] = meshSurfaceArea(mesh);

% Add fourth vertex to mesh triangles.
mesh.TRIV(:,4) = 0;

% Get ideal transformations and translations from bone transformation.
[MsIdeal,~] = Yan3DBoneTransformations(Skel1, Skel2, true);

u1 = verts(:,mesh.TRIV(:,1));
u2 = verts(:,mesh.TRIV(:,2));
u3 = verts(:,mesh.TRIV(:,3));
a = u2 - u1;
b = u3 - u2;
c = cross(a,b);
d = repmat(sum(c.^2),3,1);
u4 = (u1+u2+u3)./3 + c .* sqrt(d);
Us = [reshape(u1-u4,3,1,size(u1,2)),...
    reshape(u2-u4,3,1,size(u1,2)),...
    reshape(u3-u4,3,1,size(u1,2))];
Us = reshape(cell2mat(arrayfun(@(M)inverseRow(Us,M),1:size(Us,3),'UniformOutput',false)), size(Us,1), size(Us,2), size(Us,3));

mesh.X = [mesh.X;u4(1,:)'];
mesh.Y = [mesh.Y;u4(2,:)'];
mesh.Z = [mesh.Z;u4(3,:)'];
mesh.TRIV(:,4) = nPnts+1:numel(mesh.X);
nPnts = numel(mesh.X);

save tmp.mat

% Initialise K, dx, dy, and dz.
K = sparse(nPnts,nPnts);
% K = zeros(nPnts,nPnts);
dx = zeros(nPnts,1);
dy = zeros(nPnts,1);
dz = zeros(nPnts,1);


idx1 = mesh.TRIV(:,1);
idx2 = mesh.TRIV(:,2);
idx3 = mesh.TRIV(:,3);
idx4 = mesh.TRIV(:,4);

a = Us(1,1,:);  a=a(:);
b = Us(1,2,:);  b=b(:);
c = Us(1,3,:);  c=c(:);
d = Us(2,1,:);  d=d(:);
e = Us(2,2,:);  e=e(:);
f = Us(2,3,:);  f=f(:);
g = Us(3,1,:);  g=g(:);
h = Us(3,2,:);  h=h(:);
i = Us(3,3,:);  i=i(:);

Ma = MsIdeal(1,1,Skel1.TBA(:));  Ma=Ma(:);
Mb = MsIdeal(1,2,Skel1.TBA(:));  Mb=Mb(:);
Mc = MsIdeal(1,3,Skel1.TBA(:));  Mc=Mc(:);
Md = MsIdeal(2,1,Skel1.TBA(:));  Md=Md(:);
Me = MsIdeal(2,2,Skel1.TBA(:));  Me=Me(:);
Mf = MsIdeal(2,3,Skel1.TBA(:));  Mf=Mf(:);
Mg = MsIdeal(3,1,Skel1.TBA(:));  Mg=Mg(:);
Mh = MsIdeal(3,2,Skel1.TBA(:));  Mh=Mh(:);
Mi = MsIdeal(3,3,Skel1.TBA(:));  Mi=Mi(:);

% Precompute multiplications.
aa = a.^2;
bb = b.^2;
cc = c.^2;
dd = d.^2;
ee = e.^2;
ff = f.^2;
gg = g.^2;
hh = h.^2;
ii = i.^2;

ad = a .* d;
ag = a .* g;
be = b .* e;
bh = b .* h;
cf = c .* f;
ci = c .* i;
dg = d .* g;
eh = e .* h;
fi = f .* i;

aMa = a .* Ma;
aMd = a .* Md;
aMg = a .* Mg;

bMb = b .* Mb;
bMe = b .* Me;
bMh = b .* Mh;

cMc = c .* Mc;
cMf = c .* Mf;
cMi = c .* Mi;

dMa = d .* Ma;
dMd = d .* Md;
dMg = d .* Mg;

eMb = e .* Mb;
eMe = e .* Me;
eMh = e .* Mh;

fMc = f .* Mc;
fMf = f .* Mf;
fMi = f .* Mi;

gMa = g .* Ma;
gMd = g .* Md;
gMg = g .* Mg;

hMb = h .* Mb;
hMe = h .* Me;
hMh = h .* Mh;

iMc = i .* Mc;
iMf = i .* Mf;
iMi = i .* Mi;

% Using face area.
K = K + sparse(idx1,idx1,A.*(aa + bb + cc),nPnts,nPnts);
K = K + sparse(idx1,idx2,A.*(ad + be + cf),nPnts,nPnts);
K = K + sparse(idx1,idx3,A.*(ag + bh + ci),nPnts,nPnts);
K = K + sparse(idx1,idx4,A.*(- aa - bb - cc ...
    - ad - be - cf - ag - bh - ci),nPnts,nPnts);

dx = dx + sparse(idx1,1,A.*(- aMa - bMb - cMc),nPnts,1);
dy = dy + sparse(idx1,1,A.*(- aMd - bMe - cMf),nPnts,1);
dz = dz + sparse(idx1,1,A.*(- aMg - bMh - cMi),nPnts,1);

K = K + sparse(idx2,idx1,A.*(ad + be + cf),nPnts,nPnts);
K = K + sparse(idx2,idx2,A.*(dd + ee + ff),nPnts,nPnts);
K = K + sparse(idx2,idx3,A.*(dg + eh + fi),nPnts,nPnts);
K = K + sparse(idx2,idx4,A.*(- ad - dd - be ...
    - ee - cf - ff - dg - eh - fi),nPnts,nPnts);

dx = dx + sparse(idx2,1,A.*(- dMa - eMb - fMc),nPnts,1);
dy = dy + sparse(idx2,1,A.*(- dMd - eMe - fMf),nPnts,1);
dz = dz + sparse(idx2,1,A.*(- dMg - eMh - fMi),nPnts,1);

K = K + sparse(idx3,idx1,A.*(ag + bh + ci),nPnts,nPnts);
K = K + sparse(idx3,idx2,A.*(dg + eh + fi),nPnts,nPnts);
K = K + sparse(idx3,idx3,A.*(gg + hh + ii),nPnts,nPnts);
K = K + sparse(idx3,idx4,A.*(- ag - dg - gg ...
    - bh - eh - hh - ci - fi - ii),nPnts,nPnts);

dx = dx + sparse(idx3,1,A.*(- gMa - hMb - iMc),nPnts,1);
dy = dy + sparse(idx3,1,A.*(- gMd - hMe - iMf),nPnts,1);
dz = dz + sparse(idx3,1,A.*(- gMg - hMh - iMi),nPnts,1);

K = K + sparse(idx4,idx1,A.*(- aa - bb - cc ...
    - ad - be - cf - ag - bh - ci),nPnts,nPnts);
K = K + sparse(idx4,idx2,A.*(- dd - ee - ff ...
    - ad - be - cf - dg - eh - fi),nPnts,nPnts);
K = K + sparse(idx4,idx3,A.*(- gg - hh - ii ...
    - ag - dg - bh - eh - ci - fi),nPnts,nPnts);
K = K + sparse(idx4,idx4,A.*(aa + bb + cc ...
    + 2.*ad + dd + 2.*be + ee + 2.*cf + ff ...
    + 2.*ag + 2.*dg + gg + 2.*bh + 2.*eh + hh + 2.*ci + 2.*fi + ii),nPnts,nPnts);

dx = dx + sparse(idx4,1,A.*(aMa + dMa + gMa + bMb ...
    + eMb + hMb + cMc + fMc + iMc),nPnts,1);
dy = dy + sparse(idx4,1,A.*(aMd + dMd + gMd + bMe ...
    + eMe + hMe + cMf + fMf + iMf),nPnts,1);
dz = dz + sparse(idx4,1,A.*(aMg + dMg + gMg + bMh ...
    + eMh + hMh + cMi + fMi + iMi),nPnts,1);



% Solve linear equations.
% K = sparse(K);
X = K\dx;
Y = K\dy;
Z = K\dz;

mesh.X = X(1:nOrigPnts);
mesh.Y = Y(1:nOrigPnts);
mesh.Z = Z(1:nOrigPnts);

mesh.TRIV = mesh.TRIV(:,1:3);

return;

function M2 = inverseRow(M,r)
M2 = M(:,:,r)^-1;
return