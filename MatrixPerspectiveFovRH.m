function [matrix] = MatrixPerspectiveFovRH(fov, z_near, z_far)
    %MATRIXPERSPECTIVEFOVRH Return a matrix for the perspective view transform
    %   fov: Field-of-View in degrees
    %   z_near: Position of the near plane
    %   z_far: Position of the far place
    arguments
        fov double;
        z_near double {mustBeGreaterThan(z_near, 0)};
        z_far double {mustBeGreaterThan(z_far, z_near)}
    end

    fovy = deg2rad(fov);

    yScale = cot(fovy / 2);
    xScale = yScale;    % aspect ratio = 1 
    z_diff = z_near - z_far;

    % construct the matrix
    matrix = diag([xScale yScale]);
    matrix(3:4, 3:4) = [z_far / z_diff, -1;
                       z_near * z_far / z_diff, 0];

end
