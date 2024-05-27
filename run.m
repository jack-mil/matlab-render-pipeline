%% DEFINE PARAMETERS
LIGHT_LOC = [0,10,0];       % 10 units above (Y)
CAM_LOC = [0,1,-2];      
CAM_TARGET = [0, 0, 0];   

OBJ_LOC = [0, 0, 0];  % Location of the object
OBJ_ROT = deg2rad([0, 20, 0]);
OBJ_SCALE = 1;
OBJ_RGB = [255, 0, 0];

FOV = 90;
Z_NEAR = 0.1;
Z_FAR = 1;

% Using Right Handed Convetions

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

TL = stlread('suzanne.stl');
points = resize(TL.Points', 4, FillValue=1)';

% Do the world transformation
points_T = points * M_Scale * M_Translate * M_Rot;

figure(1);
trimesh(TL.ConnectivityList,points_T(:,1),points_T(:,2),points_T(:,3));

%% VIEW TRANSFORMATION
% -------------------

M_View = MatrixLookAtRH(CAM_LOC, CAM_TARGET);
points_view = points_T * M_View;

figure(2)
trimesh(TL.ConnectivityList,points_view(:,1),points_view(:,2),points_view(:,3));

%% PROJECTION TRANSFORMATION
% -------------------------

M_Proj = MatrixPerspectiveFovRH(FOV, Z_NEAR, Z_FAR);
points_proj = points_view * M_Proj;
points_proj = points_proj ./ points_proj(:,4);      % Perform the perspective divide operation

sorted_vertices = sortrows(points_proj,4);
sorted_faces = TL.ConnectivityList
figure(3)
triplot(TL.ConnectivityList,points_proj(:,1),points_proj(:,2));
axis([-1 1 -1 1]);
axis square;


