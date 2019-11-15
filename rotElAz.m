% rotElAz
%
% Rotates a point cloud around (0,0,0) by a given elevation and azimuth (in
% that order). z is the vertical axis.
%
% Input variables:
% pts = 3 x N matrix of x,y,z coordinates
% azel = [azimuth elevation] of the rotation
%
% Output variables:
% pts = 3 x N matrix of rotated x,y,z coordinates

function pts = rotElAz(pts,azel)

Ry = [cos(azel(2)) 0 sin(azel(2)); 0 1 0; -sin(azel(2)) 0 cos(azel(2))];
pts = Ry*pts;
Rz = [cos(-azel(1)) -sin(-azel(1)) 0; sin(-azel(1)) cos(-azel(1)) 0; 0 0 1];
pts = Rz*pts;

end