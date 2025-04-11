function Rz = rot_z(alpha)
%ROTZ Rotation matrix around z-axis
%   Rz = ROTZ(alpha) returns the 3x3 rotation matrix around z-axis.
    Rz = [cos(alpha) -sin(alpha) 0;
        sin(alpha) cos(alpha) 0;
        0 0 1];
end

