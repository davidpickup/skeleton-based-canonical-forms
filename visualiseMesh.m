function visualiseMesh(mesh, options)
% visualiseMesh(mesh, colour)
% Visualises a 3D mesh.
% Variables:
% mesh - input mesh to be visualised.
% options(optional) -   vertexColour = colour of vertices,
%                       faceColour = colour of faces,
%                       alpha = opacity.
%
% David Pickup 2015

alpha = 1;
CMap = 'jet';
if nargin == 2
    if isfield(options,'CMap')
        CMap = options.CMap;
    end
    if isfield(options,'alpha')
        alpha = options.alpha;
    end
    if isfield(options,'faceColour')
        h = patch('Faces',mesh.TRIV,'Vertices',[mesh.X,mesh.Y,mesh.Z],...
            'FaceVertexCData',options.faceColour,'FaceColor','flat',...
            'EdgeAlpha',0,'FaceAlpha',alpha);
        colormap(CMap);
    elseif isfield(options,'vertexColour')
        h = trisurf(mesh.TRIV,mesh.X,mesh.Y,mesh.Z,options.vertexColour,'FaceAlpha',alpha);
        shading flat
        colormap(CMap);
    else
        h = trisurf(mesh.TRIV,mesh.X,mesh.Y,mesh.Z,'FaceAlpha',alpha);
        shading interp, colormap([0.8, 0.559, 0.302]), lighting phong, camlight head
    end
else
    h = trisurf(mesh.TRIV,mesh.X,mesh.Y,mesh.Z,'FaceAlpha',alpha);
    shading interp, colormap([0.8, 0.559, 0.302]), lighting phong, camlight head
end

%colormap([1 1 1]*0.9), lighting phong, camlight head
set(gcf,'Renderer','OpenGL');
set(gcf,'color','white')
axis image; axis off
title('Mesh');
drawnow

return;