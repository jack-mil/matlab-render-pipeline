%% DEFINE PARAMETERS
LIGHT_LOC = [0,10,-5];
CAM_LOC = [0, 1, -3];      
CAM_TARGET = [0, 0, 0];   

OBJ_LOC = [0, 0, 0];  % Location of the object
OBJ_ROT = deg2rad([0, 0, 0]);
OBJ_SCALE = 0.5;

% OBJ_RGB = [0.3010 0.7450 0.9330];  % Blue
OBJ_RGB = [0.8500 0.3250 0.0980];  % Orange
% OBJ_RGB = [0.4660 0.6740 0.1880];  % Green
% OBJ_RGB = [0.6350 0.0780 0.1840];  % Red

FOV = 70;
Z_NEAR = 0.1;
Z_FAR = 1000;

%% LOAD MODEL

TL = stlread('cylinder_hi_poly.stl');
points = resize(TL.Points', 4, FillValue=1)';   % Add the extra w = 1 dimention to make transformation easier
tris = TL.ConnectivityList;

%% WORLD TRANSFORMATION
% --------------------
M_Scale = OBJ_SCALE * eye(3,3);
M_Scale(4,4) = 1;

M_Translate = eye(4,4);
M_Translate(4,1:3) = OBJ_LOC;    % Scale not applied to translation

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
colors = zeros(length(points_T),3);
% normals = zeros(length(tris),3);
% lights = zeros(length(tris),3);
% centers = zeros(length(tris),3);

TL = triangulation(tris, points_T(:,1:3));
V_norms = vertexNormal(TL);     % norm of each vertex (rows)
L = LIGHT_LOC - points_T(:,1:3); % light vector for every vertex
L = L ./ vecnorm(L, 2, 2);       % normalize each vector (euclidean, dim=2)
tmp = dot(L, V_norms, 2);   % dot product of the row vectors (dim=2)
Cdif = OBJ_RGB .* max(0,tmp);   % This is the vertex color

% for i = 1:size(tris, 1)
%     % Get the vertex IDs for the ith facet
%     vertex_ids = tris(i, :);
% 
%     % 3x3 matrix, rows are the x,y,z,w of each vertex
%     ver = points(vertex_ids, :);
%     ver = resize(ver,[3,3]); % get rid of the ws (1)
% 
%     % surface normal for this tri
%     N = cross(ver(2,:) - ver(1,:), ver(3,:) - ver(1,:));
%     N = N / norm(N);
%     normals(i, :) = N;
% 
%     center = mean(ver);
%     centers(i, :) = center;
%     L = LIGHT_LOC - center;
%     L = L / norm(L);
%     lights(i,:) = L;
% 
%     % V = CAM_LOC - center(1:3)
%     % V = V / norm(V);
%     tmp = dot(N, L);
%     Cdif = OBJ_RGB*max(0, tmp);
%     colors(i, :) = Cdif;
% end

figure(2)

trisurf(TL, FaceColor=[0.8,0.8,1])
% axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
axis square;
hold on
% centers = num2cell(centers, 1);
% lights = num2cell(lights, 1);
% normals = num2cell(normals, 1);


V_norms = num2cell(V_norms, 1);
quiver3(points_T(:,1),points_T(:,2),points_T(:, 3), ...
    V_norms{:}, color='g');
% quiver3(centers{:}, ...
%     normals{:}, Color='r')
% quiver3(centers{:}, ...
%     lights{:}, Color='b')
hold off

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
points_proj = points_proj ./ points_proj(:, 4);      % Perform the perspective divide operation

%% Z SORTING
% ----------
% do some fancy z-sorting for painter's algo
avg_z = mean(points_proj(tris, 3), 2);
zs = points_proj(tris(:), 3); % get the z coord for each vertex for every triangle
zs = reshape(zs, size(tris));    % and group each z position by triangle

avg_z = mean(zs, 2); % column vector with z of each triangle

% combine z-location, vertices for a tri, and color into a Nx7 matrix.
% This keeps all the information about a particular tri together when sorting
% aug = [avg_z tris colors];
aug = [avg_z, tris];        % no need to associate colors with facets when using vertex shading

sorted_tris = sortrows(aug, 1, "descend");   % sort based on first column


%% RESTERIZATION
% --------------
figure(3)
clf(3)

% patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", sorted_tris(:,5:7), "FaceColor", "flat", "EdgeColor", "none");
patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Cdif, "FaceColor", "interp", "EdgeColor", "none");

axis([-1 1 -1 1]);
axis square;


