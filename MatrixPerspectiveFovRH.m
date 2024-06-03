function [matrix] = MatrixPerspectiveFovRH(fov, z_near, z_far, wh_ratio)
    %MATRIXPERSPECTIVEFOVRH Return a matrix for the perspective view transform
    %   fov: Field-of-View in degrees
    %   z_near: Position of the near plane
    %   z_far: Position of the far place
    %   wh_ratio: Vector with relative size of axis [w h]
    arguments
        fov double;
        z_near double {mustBeGreaterThan(z_near, 0)};
        z_far double {mustBeGreaterThan(z_far, z_near)}
        wh_ratio(1, 2) double = [1, 1];
    end

    fovy = deg2rad(fov);

    yScale = cot(fovy / 2);
    xScale = yScale / (wh_ratio(1) / wh_ratio(2));
    z_diff = z_near - z_far;

    % construct the matrix
    matrix = diag([xScale yScale]);
    matrix(3:4, 3:4) = [z_far / z_diff, -1;
                       z_near * z_far / z_diff, 0];

end
