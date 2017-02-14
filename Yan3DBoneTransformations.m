function [Ms, Ts] = Yan3DBoneTransformations(Skel1, Skel2, askForArray)
% [Ms, Ts] = Yan3DBoneTransformations(Skel1, Skel2)
% Calculates the shape transformation matrix and translation vector from
% each bone in Skel1 to Skel2, using the method by Yan et al. TVCG 2008.
% Variables:
% Ms - shape transformation matrices.
% Ts - translation matrices.
% Skel1 - structure of first skeleton.
% Skel2 - structure of second skeleton.
%
% David Pickup 2014

if nargin == 2
    askForArray = false;
end

% Get number of bones.
nBones = size(Skel1.E,1);

% Initialise list of transformation matrices and translation vectors.
if askForArray
    Ms = zeros(3,3,nBones);
else
    Ms = cell(nBones,1);
end
Ts = cell(nBones,1);

% Iterate through each bone.
for i = 1:nBones
    % Get verticex indices for current bone in each skeleton.
    A1 = Skel1.E(i,1);
    B1 = Skel1.E(i,2);
    A2 = Skel2.E(i,1);
    B2 = Skel2.E(i,2);
    
    % Calculate rotation translations for bone.
    TR1 = [-Skel1.X(A1); -Skel1.Y(A1); -Skel1.Z(A1)];
    TR2 = [Skel2.X(A2); Skel2.Y(A2); Skel2.Z(A2)];
    
    % Calculate rotation for bone.
    V1 = [Skel1.X(B1)-Skel1.X(A1); Skel1.Y(B1)-Skel1.Y(A1); Skel1.Z(B1)-Skel1.Z(A1)];
    V1 = V1 ./ norm(V1);
    V2 = [Skel2.X(B2)-Skel2.X(A2); Skel2.Y(B2)-Skel2.Y(A2); Skel2.Z(B2)-Skel2.Z(A2)];
    V2 = V2 ./ norm(V2);
    N = cross(V1,V2);
    a = N(1);
    b = N(2);
    c = N(3);
    theta = acos(dot(V1,V2));
    if ~isreal(theta)
        theta = 0;
    end
    mu = sin(theta);
    nu = cos(theta);
    lamda = (1-nu);
    rho_ab = a^2 + b^2;
    rho_bc = b^2 + c^2;
    rho_ac = a^2 + c^2;
    
    
%     % Rotation matrix presented in paper.
%     R = [a^2+rho_bc*nu    a*b*lamda+c*mu    a*c*lamda-b*mu;
%         a*b*lamda-c*mu    b^2+rho_ac*nu    b*c*lamda-a*b*mu;
%         a*c*lamda+b*mu    b*c*lamda-a*b*mu    c^2+rho_ab*nu];

    % Rotation matrix presented on Wikipedia.
    R = [nu+a^2*(lamda)    a*b*lamda-c*mu    a*c*lamda+b*mu;
        b*a*lamda+c*mu    nu+b^2*lamda    b*c*lamda-a*mu;
        c*a*lamda-b*mu    c*b*lamda+a*mu    nu+c^2*lamda];
    R = vrrotvec2mat([N;theta]);
    
    % Store current transformation and translation.
    if askForArray
        Ms(:,:,i) = R;
    else
        Ms{i} = R;
    end
    Ts{i} = R*TR1 + TR2;
end

return;