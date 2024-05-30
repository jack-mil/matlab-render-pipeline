%% DEFINE PARAMETERS
LIGHT_LOC = [-3,10,-5];
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

TL = stlread('teapot.stl');
points = resize(TL.Points', 4, FillValue=1)';   % Add the extra w = 1 dimention to make transformation easier
tris = TL.ConnectivityList;

%% WORLD TRANSFORMATION
% --------------------
M_Scale = OBJ_SCALE * eye(3,3);
M_Scale(4,4) = 1;

M_Translate = eye(4,4);
M_Translate(4,1:3) = OBJ_LOC;    % Scale not applied to translation

a = OBJ_ROT(1);

rotx = eye(4,4);
rotx(2:3,2:3) = [cos(a) -sin(a);
                 sin(a)  cos(a)];

a = OBJ_ROT(2);
roty = [cos(a) 0 sin(a);
    0 1 0;
    -sin(a) 0 cos(a)];
roty(4,4) = 1;

a = OBJ_ROT(3);
rotz = eye(4,4);
rotz(1:2,1:2) = [cos(a) -sin(a);
                 sin(a) cos(a)];

M_Rot = rotx * roty * rotz;

% Do the world transformation
points_T = points * M_Scale * M_Translate * M_Rot;


% figure(1);
% trimesh(TL.ConnectivityList,points_T(:, 1),points_T(:, 2),points_T(:, 3));

%% LIGHTING CALCULATIONS
% ----------------------

% colum vector to hold the calculated colors
colors = zeros(length(tris),3);
normals = zeros(length(tris),3);
lights = zeros(length(tris),3);
centers = zeros(length(tris),3);
for i = 1:size(tris, 1)
    % Get the vertex IDs for the ith facet
    vertex_ids = tris(i, :);

    % 3x3 matrix, rows are the x,y,z,w of each vertex
    ver = points(vertex_ids, :);
    ver = resize(ver,[3,3]); % get rid of the ws (1)
    
    % surface normal for this tri
    N = cross(ver(2,:) - ver(1,:), ver(3,:) - ver(1,:));
    N = N / norm(N);
    normals(i, :) = N;

    center = mean(ver);
    centers(i, :) = center;
    L = LIGHT_LOC - center;
    L = L / norm(L);
    lights(i,:) = L;

    % V = CAM_LOC - center(1:3)
    % V = V / norm(V);
    tmp = dot(N, L);
    Cdif = OBJ_RGB*max(0, tmp);
    colors(i, :) = Cdif;

end

figure(2)
trisurf(TL,FaceColor=[0.8,0.8,1])
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
axis square;
hold on
centers = num2cell(centers, 1);
lights = num2cell(lights, 1);
normals = num2cell(normals, 1);
quiver3(centers{:}, ...
    normals{:}, Color='r')
quiver3(centers{:}, ...
    lights{:}, Color='b')
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

% do some fancy z-sorting
avg_z = mean(points_proj(tris, 3), 2);
zs = points_proj(tris(:), 3); % get the z coord for each vertex for every triangle
zs = reshape(zs, size(tris));    % and group each z position by triangle

avg_z = mean(zs,2); % column vector with z of each triangle

aug = [avg_z tris colors];

sorted_tris = sortrows(aug, 1, "descend");   % sort based on first column

figure(3)
clf(3)

patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", sorted_tris(:,5:7), "FaceColor", "flat", "EdgeColor", "none");

axis([-1 1 -1 1]);
axis square;


