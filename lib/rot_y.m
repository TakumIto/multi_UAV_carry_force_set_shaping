function Ry = rot_y(beta)
%ROTY Rotation matrix around y-axis
%   Ry = ROTY(beta) returns the 3x3 rotation matrix around y-axis.
    Ry = [cos(beta) 0 sin(beta);
        0 1 0;
        -sin(beta) 0 cos(beta)];
end

