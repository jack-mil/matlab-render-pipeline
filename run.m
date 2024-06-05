%% DEFINE PARAMETERS
% ------------------
model_file = 'suzanne_hip.stl';
% Object transform
OBJ_LOC = [0, 0, 0];
OBJ_ROT = deg2rad([-90, 0, 180]);
OBJ_SCALE = 1.333;

% Object material
OBJ_RGB = [0.8,0.8,0.8];    % White
% OBJ_RGB = [0.3010 0.7450 0.9330]; % Blue
% OBJ_RGB = [1, 0.8393, 0];     % Gold
% OBJ_RGB = [0.8500 0.3250 0.0980];  % Orange
% OBJ_RGB = [0.4660 0.6740 0.1880];  % Green
% OBJ_RGB = [0.6350 0.0780 0.1840]; % Red
OBJ_Ks = OBJ_RGB; % Color of the specular highlights
% OBJ_Ks = [1,1,1];
OBJ_spec = 150; % Shininess constant

% Scene lighting
LIGHT1_LOC = [5, 3, 5]; % Key light
LIGHT1_RGB = [1, 1, 1];

LIGHT2_LOC = [-5, 2, 4]; % Fill light
LIGHT2_RGB = [0.8, 0.8, 0.8];

LIGHT3_LOC = [-5, 3, -5]; % Backlight
% LIGHT3_RGB = [0.3, 0.3, 0.3];
LIGHT3_RGB = [0.6, 0, 0.4];     % Weird reddish underglow

LIGHT4_LOC = [5,3,-5]
LIGHT4_RGB = [0,0.6,0.4];

% Add all the lights and colors into multidimensional arrays
Ls = cat(3, LIGHT1_LOC, LIGHT2_LOC, LIGHT3_LOC, LIGHT4_LOC);
L_RGBs = cat(3, LIGHT1_RGB, LIGHT2_RGB, LIGHT3_RGB, LIGHT4_RGB);

Camb = [0, 0, 0]; % No ambient light

% Camera settings
orbit = @(a, d, y) [d * cos(deg2rad(270 + a)), y, d * sin(deg2rad(270 + a))];
% CAM_LOC = orbit(0, 5, 0);
CAM_LOC = [0, 0, 5];
CAM_TARGET = [0, 0, 0];
FOV = 39.6; % 50mm
WH_RATIO = [4, 3]; % Aspect ratio
Z_NEAR = 0.1;
Z_FAR = 1000;

%% LOAD MODEL
% -----------
TL = stlread(model_file);
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
points_T = points * M_Rot * M_Scale * M_Translate;

%% LIGHTING CALCULATIONS (in world space)
% ----------------------

% Gouraud shading with Phong reflection model from
% https://en.wikipedia.org/wiki/Phong_reflection_model
% https://en.wikipedia.org/wiki/Gouraud_shading

% Useful function handles
normr = @(M) M ./ vecnorm(M, 2, 2); % Euclidean normalize every row
dotr = @(T, M) sum(T .* M, 2); % Take a dot product of rows in M and T, across pages of T

% Reconstruct the tri object for vertexNormal function
N = vertexNormal(triangulation(tris, points_T(:, 1:3))); % Normal of each vertex (rows)

L = normr(Ls - points_T(:, 1:3)); % Unit light vector of each vertex
V = normr(CAM_LOC - points_T(:, 1:3)); % Unit view vector of each vertex

LN_dot = dotr(L, N);

Cdiff = L_RGBs .* OBJ_RGB .* max(0, LN_dot); % Calculate diffuse color

R = 2 .* LN_dot .* N - L; % Reflection vector

% Calculate specular highlights
% Specular only present when diffuse is non-zero
Cspec = L_RGBs .* OBJ_Ks .* max(0, dotr(R, V) .^ OBJ_spec) .* (Cdiff > 0);

% Add the effect of all lights together (3-dimension array)
Ctot = Camb + sum(Cdiff, 3) + sum(Cspec, 3);

%% VIEW TRANSFORMATION (RH)
% -------------------

points_view = points_T * MatrixLookAtRH(CAM_LOC, CAM_TARGET);

%% PROJECTION TRANSFORMATION (RH)
% -------------------------

% Perform the perspective divide operation
points_proj = points_view * MatrixPerspectiveFovRH(FOV, Z_NEAR, Z_FAR, WH_RATIO);
points_proj = points_proj ./ points_proj(:, 4);
% Optional ortho projection
% points_proj = points_view * MatrixOrthoRH(WH_RATIO(1), WH_RATIO(2), Z_NEAR, Z_FAR);

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
figure(1);
patch('Faces', sorted_tris(:, 2:4), "Vertices", points_proj(:, 1:2), "FaceVertexCData", Ctot, "FaceColor", "interp", "EdgeColor", "none");

% Configure figure and plot dimensions
set(gca, "Color", [0.4, 0.4, 0.4]);
set(gcf, "Color", [0.8, 0.8, 0.8]);
axis([-1 1 -1 1]); pbaspect([WH_RATIO, 1])
