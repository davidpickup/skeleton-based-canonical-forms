function TBA = triangleBoneAssignment(Skel, mesh)
% TBA = triangleBoneAssignment(Skel, mesh)
% Assigns each triangle of the mesh a control bone.
% Variables:
% TBA - bone assignment of each triangle.
% Skel - skeleton structure.
% mesh - mesh structure.
%
% David Pickup 2014

%% New Way

nTris = size(mesh.TRIV,1);
nPnts = numel(mesh.X);
nJnts = numel(Skel.X);
nBns = size(Skel.E,1);

% Compute vertex-joint assignment.
VJA = zeros(nPnts,1);
for i = 1:nJnts
    VJA(Skel.H{i}) = i;
end

% Compute triangle-bone assignment.
TBA = zeros(nTris,1);
for i = 1:nJnts
    % Get bones connected to current joint.
    [B,~] = find(Skel.E==i);
    
    % Get neighbouring joint for each bone.
    N = zeros(numel(B),1);
    for j = 1:numel(B)
        if (Skel.E(B(j),1) ~= i)
            N(j) = Skel.E(B(j),1);
        else
            N(j) = Skel.E(B(j),2);
        end
    end
    
    % Get all triangles with a vertex assigned to the joint.
    [I,~] = find(ismember(mesh.TRIV,Skel.H{i}));
    I = unique(I);
    
    % If joint is an end joint.
    if numel(B) == 1
        % Assign triangles to the bone connected to the joint.
        TBA(I) = B;
    else
        % Iterate through all neighbouring joints.
        for j = 1:numel(N)
            % Compute vector in direction of neighbour.
            V1 = [  Skel.X(N(j)) - Skel.X(i);...
                    Skel.Y(N(j)) - Skel.Y(i);...
                    Skel.Z(N(j)) - Skel.Z(i)];
            V1 = V1 ./ norm(V1);
            % Initialise list of plane normals.
            Ps = zeros(3,numel(N)-1);
            % Compute bisection plane normals.
            count = 1;
            for k = 1:numel(N)
                if k~=j
                    % Compute vector in direction of second neighbour.
                    V2 = [  Skel.X(N(k)) - Skel.X(i);...
                            Skel.Y(N(k)) - Skel.Y(i);...
                            Skel.Z(N(k)) - Skel.Z(i)];
                    V2 = V2 ./ norm(V2);
                    bisection = (V1+V2) ./ norm(V1+V2);
                    perpendicular = cross(V1,V2);
                    normal = cross(bisection,perpendicular);
                    Ps(:,count) = normal ./ norm(normal);
                    count = count + 1;
                end
            end
            
            % Iterate through all triangles connected to joint by a vertex.
            for k = 1:numel(I)
                % If the triangle is already assigned, then skip.
                if TBA(I(k)) > 0
                    continue;
                end
                
                % Calculate vector to each triangle vertex.
                Vs = zeros(3,3);
                
                Vs(:,1) = [ mesh.X(mesh.TRIV(I(k),1)) - Skel.X(i);...
                            mesh.Y(mesh.TRIV(I(k),1)) - Skel.Y(i);...
                            mesh.Z(mesh.TRIV(I(k),1)) - Skel.Z(i)];
                Vs(:,1) = Vs(:,1) ./ norm(Vs(:,1));
                Vs(:,2) = [ mesh.X(mesh.TRIV(I(k),2)) - Skel.X(i);...
                            mesh.Y(mesh.TRIV(I(k),2)) - Skel.Y(i);...
                            mesh.Z(mesh.TRIV(I(k),2)) - Skel.Z(i)];
                Vs(:,2) = Vs(:,2) ./ norm(Vs(:,2));
                Vs(:,3) = [ mesh.X(mesh.TRIV(I(k),3)) - Skel.X(i);...
                            mesh.Y(mesh.TRIV(I(k),3)) - Skel.Y(i);...
                            mesh.Z(mesh.TRIV(I(k),3)) - Skel.Z(i)];
                Vs(:,3) = Vs(:,3) ./ norm(Vs(:,3));
                
                tmp = repmat(Vs,1,1,size(Ps,2));
                tmp2 = repmat(Ps,1,1,3);
                tmp3 = permute(tmp2,[1 3 2]);
                tmp4 = tmp.*tmp3;
                tmp5 = sum(tmp4,1);
                tmp5 = reshape(tmp5,[3,size(Ps,2)])';
                tmp6 = tmp5>0;
                flags = sum(tmp6,1)==size(Ps,2);
                if sum(flags) >= 2
                    TBA(I(k)) = B(j);
                end
            end
        end
    end
end

% Assign any unassigned triangles to any bone assigned to any of its
% neighbours.
Idx = find(TBA == 0);
while ~isempty(Idx)
    for i = 1:numel(Idx)
        [I,~] = find(ismember(mesh.TRIV, mesh.TRIV(Idx(i),:)));
        I = unique(I);
        X = find(TBA(I) > 0);
        if ~isempty(X)
            TBA(Idx(i)) = TBA(I(X(randi(numel(X)))));
        end
    end
    Idx = find(TBA == 0);
end



return;