function [matrix] = MatrixOrthoRH(w, h, z_near, z_far)
    %MATRIXPERSPECTIVEFOVRH Return a matrix for the orthographic projection transform
    %   w: Width of the view volume
    %   h: Height of the view volume
    %   z_near: Position of the near plane
    %   z_far: Position of the far place
    arguments
        w double {mustBeGreaterThan(w, 0)};
        h double {mustBeGreaterThan(h, 0)};
        z_near double {mustBeGreaterThan(z_near, 0)};
        z_far double {mustBeGreaterThan(z_far, z_near)}
    end

    z_diff = z_near - z_far;
    % Construct the matrix
    matrix = diag([2/w 2/h]);
    matrix(3:4, 3:4) = [1 / z_diff, 0;
                       z_near / z_diff, 1];

end
