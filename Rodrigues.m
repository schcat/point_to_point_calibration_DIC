function [rotation_matrix] = Rodrigues( rotation_vector )
%   Convert a rotation vector to a rotation matrix.
    rotation_vector = rotation_vector';
    theta = norm(rotation_vector);
    rotation_vector = rotation_vector./theta;
    I = eye(3);
    tmp_matrix = [0 -rotation_vector(3) rotation_vector(2); ...
        rotation_vector(3) 0 -rotation_vector(1); ...
        -rotation_vector(2) rotation_vector(1) 0];
    rotation_matrix = cos(theta) * I + (1 - cos(theta)) * ...
    (rotation_vector * rotation_vector') + sin(theta) * tmp_matrix;
end