%% DEFINE PARAMETERS
% ------------------
% Scene lighting
LIGHT_LOC = [-10, 3, -4.5]; % Simple directional light (sun)
% LIGHT_RGB = [1, 1, 1];      % White Light
% LIGHT_RGB = [1, 0.75, 0.8]; % Soft Pink Light
% LIGHT_RGB = [1, 0, 0];      % Red light
% LIGHT_RGB = [1, 1, 0];      % Yellow light
% Camb = [0.2, 0.2, 0.2];

% Camera settings
CAM_LOC = [0, 2, -5];
CAM_TARGET = [0, 0, 0];
FOV = 70;
Z_NEAR = 0.1;
Z_FAR = 1000;

% Object transform
OBJ_LOC = [0, 0, 0];
OBJ_ROT = deg2rad([0, 0, 0]);
OBJ_SCALE = 1;

% OBJ_RGB = [1,1,1];     % Grey / silver
OBJ_RGB = [0.3010 0.7450 0.9330]; % Blue
% OBJ_RGB = [0.8500 0.3250 0.0980];  % Orange
% OBJ_RGB = [0.4660 0.6740 0.1880];  % Green
% OBJ_RGB = [0.6350 0.0780 0.1840];  % Red
OBJ_Ks = [1, 1, 1]; % color of the specular light (gray
OBJ_spec = 50; % shininess
Camb = OBJ_RGB .* 0.4;

% Useful function handles
normr = @(M) M ./ vecnorm(M, 2, 2); % euclidean normalize every row

%% LOAD MODEL
% -----------

TL = stlread('sphere_hip.stl');
points = resize(TL.Points', 4, FillValue = 1)'; % Add the extra w = 1 dimenstion to make transformation easier
tris = TL.ConnectivityList;

%% WORLD TRANSFORMATION
% --------------------
M_Scale = OBJ_SCALE * eye(3, 3);
M_Scale(4, 4) = 1;

M_Translate = eye(4, 4);
M_Translate(4, 1:3) = OBJ_LOC; % Scale not applied to translation

% a = OBJ_ROT(1);
% rotx = eye(4,4);
% rotx(2:3,2:3) = [cos(a) -sin(a);
%                  sin(a)  cos(a)];
rotx = @(a) [1, 0, 0, 0; 0, cos(a), -sin(a), 0; 0, sin(a), cos(a), 0; 0, 0, 0, 1];

% a = OBJ_ROT(2);
% roty = [cos(a) 0 sin(a);
%     0 1 0;
%     -sin(a) 0 cos(a)];
% roty(4,4) = 1;
roty = @(a) [cos(a), 0, sin(a), 0; 0, 1, 0, 0; -sin(a), 0, cos(a), 0; 0, 0, 0, 1];

% a = OBJ_ROT(3);
rotz = @(a) [cos(a), -sin(a), 0, 0; sin(a), cos(a), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
% rotz = eye(4,4);
% rotz(1,:) = [cos(a), -sin(a), 0, 0];
% rotz(2,:) = [sin(a), cos(a), 0, 0];
% rotz(1:2,1:2) = [cos(a) -sin(a);
% sin(a) cos(a)];

% Create the overall rotation matrix
M_Rot = rotx(OBJ_ROT(1)) * roty(OBJ_ROT(2)) * rotz(OBJ_ROT(3));

% Do the world transformation
points_T = points * M_Scale * M_Translate * M_Rot;

% figure(1);
% trimesh(TL.ConnectivityList,points_T(:, 1),points_T(:, 2),points_T(:, 3));

%% LIGHTING CALCULATIONS
% ----------------------

% colum vector to hold the calculated colors
colors = zeros(length(points_T), 3);
% normals = zeros(length(tris),3);
% lights = zeros(length(tris),3);
% centers = zeros(length(tris),3);

TL = triangulation(tris, points_T(:, 1:3));
N = vertexNormal(TL); % norm of each vertex (rows)

L = LIGHT_LOC - points_T(:, 1:3); % light vector for each vertex
V = CAM_LOC - points_T(:, 1:3); % view vector for each vertex

L = normr(L); % normalize each vector (euclidean, dim=2)
V = normr(V);

LN_dot = dot(L, N, 2); % calculate angle of light and vertex normals

Cdiff = OBJ_RGB .* max(0, LN_dot);

R = 2 .* LN_dot .* N - L; % Reflection vector
Cspec = OBJ_Ks .* max(0, dot(R, V, 2) .^ OBJ_spec) .* max(0, LN_dot);

Ctot = Camb + Cdiff + Cspec;

% figure(2)
%
% trisurf(TL, FaceColor=[0.8,0.8,1])
% % axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
% axis square;
% hold on
% % centers = num2cell(centers, 1);
% % lights = num2cell(lights, 1);
% % normals = num2cell(normals, 1);
%
%
% V_norms = num2cell(V_norms, 1);
% quiver3(points_T(:,1),points_T(:,2),points_T(:, 3), ...
%     V_norms{:}, color='g');
% quiver3(centers{:}, ...
%     normals{:}, Color='r')
% quiver3(centers{:}, ...
%     lights{:}, Color='b')
% hold off

%% VIEW TRANSFORMATION
% -------------------

M_View = MatrixLookAtRH(CAM_LOC, CAM_TARGET);
points_view = points_T * M_View;

% figure(2)
% trimesh(TL.ConnectivityList, points_view(:, 1), points_view(:, 2), points_view(:, 3));

%% PROJECTION TRANSFORMATION
% -------------------------

M_Proj = MatrixPerspectiveFovRH(FOV, Z_NEAR, Z_FAR);
points_proj = points_view * M_Proj;
points_proj = points_proj ./ points_proj(:, 4); % Perform the perspective divide operation

%% Z SORTING
% ----------
% do some fancy z-sorting for painter's algo
avg_z = mean(points_proj(tris, 3), 2);
zs = points_proj(tris(:), 3); % get the z coord for each vertex for every triangle
zs = reshape(zs, size(tris)); % and group each z position by triangle

avg_z = mean(zs, 2); % column vector with z of each triangle

% combine z-location, vertices for a tri, and color into a Nx7 matrix.
% This keeps all the information about a particular tri together when sorting
% aug = [avg_z tris colors];
aug = [avg_z, tris]; % no need to associate colors with facets when using vertex shading

sorted_tris = sortrows(aug, 1, "descend"); % sort based on first column

%% RESTERIZATION
% --------------

% f1=figure(1);
% patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceColor" , Camb, "EdgeColor", "none");
% axis([-1 1 -1 1]);
% axis square;
% f2=figure(2);
% patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Cdiff, "FaceColor", "interp", "EdgeColor", "none");
% axis([-1 1 -1 1]);
% axis square;
% f3=figure(3);
% patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Cspec, "FaceColor", "interp", "EdgeColor", "none");
% axis([-1 1 -1 1]);
% axis square;
f4 = figure(4);
patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Ctot, "FaceColor", "interp", "EdgeColor", "none");
axis([-1 1 -1 1]); axis square;
