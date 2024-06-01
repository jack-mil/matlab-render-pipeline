%% DEFINE PARAMETERS
% ------------------
% Scene lighting
LIGHT_LOC = [-10, 3, -4.5]; % Simple directional light (sun)
LIGHT_RGB = [1, 1, 1];      % White Light
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

% Object material
OBJ_RGB = [0.3010 0.7450 0.9330]; % Blue
% OBJ_RGB = [1,1,1];     % Grey / silver
% OBJ_RGB = [0.8500 0.3250 0.0980];  % Orange
% OBJ_RGB = [0.4660 0.6740 0.1880];  % Green
% OBJ_RGB = [0.6350 0.0780 0.1840];  % Red
OBJ_Ks = [1, 1, 1]; % Color of the specular highlights
OBJ_spec = 50; % Shininess constant
Camb = OBJ_RGB .* 0.4;  % Ambient light color

% Useful function handles
normr = @(M) M ./ vecnorm(M, 2, 2); % Euclidean normalize every row

%% LOAD MODEL
% -----------

TL = stlread('teapot.stl');
% Add the extra w = 1 dimension to make transformation easier.
points = resize(TL.Points', 4, FillValue = 1)';
tris = TL.ConnectivityList;

%% WORLD TRANSFORMATION
% --------------------
M_Scale = OBJ_SCALE * eye(3, 3);
M_Scale(4, 4) = 1;

M_Translate = eye(4, 4);
M_Translate(4, 1:3) = OBJ_LOC; % Scale not applied to translation

% Define 4x4 rotation operators
rotx = @(a) [1, 0, 0, 0; 0, cos(a), -sin(a), 0; 0, sin(a), cos(a), 0; 0, 0, 0, 1];

roty = @(a) [cos(a), 0, sin(a), 0; 0, 1, 0, 0; -sin(a), 0, cos(a), 0; 0, 0, 0, 1];

rotz = @(a) [cos(a), -sin(a), 0, 0; sin(a), cos(a), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];

% Create the overall rotation matrix
M_Rot = rotx(OBJ_ROT(1)) * roty(OBJ_ROT(2)) * rotz(OBJ_ROT(3));

% Do the world transformation
points_T = points * M_Scale * M_Translate * M_Rot;

%% LIGHTING CALCULATIONS (in world space)
% ----------------------

% Gouraud shading with Phong reflection model from
% https://en.wikipedia.org/wiki/Phong_reflection_model

% Reconstruct the tri object for vertexNormal function
N = vertexNormal(triangulation(tris, points_T(:, 1:3))); % Norm of each vertex (rows)

L = normr(LIGHT_LOC - points_T(:, 1:3));    % Unit light vector of each vertex
V = normr(CAM_LOC - points_T(:, 1:3));      % Unit view vector of each vertex

LN_dot = dot(L, N, 2); % angle b/w/ normal and incoming light

Cdiff = LIGHT_RGB .* OBJ_RGB .* max(0, LN_dot);  % Calculate diffuse color

R = 2 .* LN_dot .* N - L; % Reflection vector

% Calculate specular highlights
% Specular only present when diffuse is non-zero
Cspec = LIGHT_RGB .* OBJ_Ks .* max(0, dot(R, V, 2) .^ OBJ_spec) .* (Cdiff>0);

Ctot = Camb + Cdiff + Cspec;

%% VIEW TRANSFORMATION
% -------------------

points_view = points_T * MatrixLookAtRH(CAM_LOC, CAM_TARGET);

%% PROJECTION TRANSFORMATION
% -------------------------

points_proj = points_view * MatrixPerspectiveFovRH(FOV, Z_NEAR, Z_FAR);
% Perform the perspective divide operation
points_proj = points_proj ./ points_proj(:, 4);

%% Z SORTING
% ----------
% Do some fancy z-sorting for painter's algo
zs = points_proj(tris(:), 3); % Get the z pos of each vertex for every facet,
zs = reshape(zs, size(tris)); % and reshape to Nx3 group each z position by triangle

avg_z = mean(zs, 2); % Get the mean of every row

sorted_tris = sortrows([avg_z, tris], 1, "descend"); % Sort based on avg_z (col 1)

%% RASTERIZATION
% --------------

% Individual passes
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
clf
f4 = figure(4);
patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Ctot, "FaceColor", "interp", "EdgeColor", "none");
axis([-1 1 -1 1]); axis square;
