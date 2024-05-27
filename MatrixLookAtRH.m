function [matrix] = MatrixLookAtRH(loc, target, up)
%MATRIXLOOKATRH Return a matrix for the camera look at transformation
%   loc: Position of the camera
%   target: Direction camera looks in
%   up: "Up" direction of the camera. Should not be colinear with target
    arguments
        loc(1,3) double;
        target(1,3) double;
        up(1,3) double = [0,1,0];
    end
    zaxis = loc - target;         % normalize
    zaxis = zaxis / norm(zaxis);
    xaxis = cross(up, zaxis);
    xaxis = xaxis / norm(xaxis);
    yaxis = cross(zaxis, xaxis);

    matrix = [xaxis; yaxis; zaxis; zeros(1,3)];
    matrix(4, 1:4) = [-dot(xaxis, loc), -dot(yaxis, loc) -dot(zaxis,loc), 1];
end