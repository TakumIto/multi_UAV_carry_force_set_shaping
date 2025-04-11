function Rx = rot_x(gamma)
%ROTX Rotation matrix around x-axis
%   Rx = ROTX(gamma) returns the 3x3 rotation matrix around x-axis.
    Rx = [1 0 0;
        0 cos(gamma) -sin(gamma);
        0 sin(gamma) cos(gamma)];
end

