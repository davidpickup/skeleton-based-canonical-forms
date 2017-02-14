function outMesh = AuMeshContraction(inMesh)
% outMesh = AUMeshContraction(inMesh)
% Contracts a mesh's geometry using the method by Au et al. SIGGRAPH 2008.
% Variables:
% outMesh - contracted mesh.
% inMesh - input mesh structure.
%
% David Pickup 2013

MaxWl = 2048;
MaxWh = 10000;

nTris = size(inMesh.TRIV,1);
nVerts = numel(inMesh.X);

% Get mesh surface area.
[TotalArea,AreaFace] = meshSurfaceArea(inMesh);
AverageArea = TotalArea / nTris;

% Initialise weights.
Wh = ones(nVerts,1);
% Wl = ones(nVerts,1).*3;
Wl = ones(nVerts,1).*((10^(-3)) * sqrt(AverageArea));

% Initialise contracted vertices.
outMesh = inMesh;
verts = [outMesh.X,outMesh.Y,outMesh.Z];

meanRatio = 0;
for i = 1:20
    L = laplaceBeltrami(inMesh);
    
    A = [sparse(1:nVerts, 1:nVerts, Wl)*L;sparse(1:nVerts, 1:nVerts, Wh)];
    B = [zeros(size(verts));sparse(1:nVerts, 1:nVerts, Wh)*verts];
    verts = A\B;

    inMesh.X = verts(:,1);
    inMesh.Y = verts(:,2);
    inMesh.Z = verts(:,3);
    
    Wl = 3.0*Wl;
    Wl(Wl>MaxWl) = MaxWl;
    
    [TotalAreaNew,AreaFaceNew] = meshSurfaceArea(inMesh);
    
    totalRatio = TotalAreaNew/TotalArea;
    
    ratio = neighbourhoodAreaRatio(AreaFace, AreaFaceNew, inMesh.TRIV, nVerts);
    meanRatio = mean(ratio);
    if isnan(meanRatio) || meanRatio==Inf || totalRatio <= 0.01
        if isnan(meanRatio)
            'NaN'
        elseif meanRatio==Inf
            'Inf'
        else
            outMesh = inMesh;
        end
        break;
    end
    
    outMesh = inMesh;
    
    Wh = Wh .* sqrt(ratio);
    Wh(Wh>MaxWh) = MaxWh;
    
    X = find(AreaFaceNew == 0);
    if numel(X) > 0
        [int2str(numel(X)) ' degenerate faces.']
        break;
    end
    
%     visualiseMesh(inMesh);
    AreaFace = AreaFaceNew;
end
[int2str(i) ' iterations.']
return;