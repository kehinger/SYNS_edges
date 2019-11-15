% pano2view
% Makes a camera-view image from a panorama specified by panoid.
%
% Input variables:
% pano = panoramic image (assumed to be double and in color)
% a = [az el] of view (az = 0-360, el = (-90)-90, az = 0 corresponds to
%   the left edge of the pano (unless using world coordinates, see below).
%   If view is a scalar it is assumed to be az and el = 0.
% fov = [horz vert] specifies FOV of the camera (in degrees)
%
% Output variables:
% img = projected view
% coords = az,el coordinates of the pixels in the projected view, in the
%   original panorama (size is imH x imW x 2)

function [img,coords] = pano2view(pano,a,fov,outSize)

% Set elevation to 0 if not specified in input
if length(a) == 1
    a = [a 0];
end

% 3d coordinates of pano pixels
th1 = repmat(0:size(pano,2)-1,[size(pano,1) 1]);
th1 = ((2*pi)/size(pano,2))*(th1+0.5);
th2 = repmat(size(pano,1)/2:-1:-size(pano,1)/2+1,[size(pano,2) 1]);
th2 = (pi/size(pano,1))*(th2'-0.5);
[x,y,z] = sph2cart(th1,th2,1);

% Rotate in 3d so that selected view is centered at (0,0)
azView = (pi/180)*a(1);
elView = (pi/180)*a(2);
px = [x(:) y(:) z(:)]';
px = rotAzEl(px,[azView elView]);

% Convert back to spherical coordinates
[theta,phi,~] = cart2sph(px(1,:),px(2,:),px(3,:));

% Get the coordinates of pixels in the selected view (this function
% automatically chooses resolution based on the size of the panorama)
[az,el] = pixelCoords(size(pano),fov,outSize/2);

% Rotate coordinates of the pixels in this view back to their positions in
% the original panorama (these will be saved and used to plot view data
% back on sphere)
[vx,vy,vz] = sph2cart(az,el,1);
vpts = [vx(:) vy(:) vz(:)]';
vpts = rotElAz(vpts,-[azView elView]);
[thOut,phiOut,~] = cart2sph(vpts(1,:),vpts(2,:),vpts(3,:));
% Make pixels go from 0-2pi instead of +/-pi (for compatibility with other code)
k = find(thOut < 0);
thOut(k) = thOut(k) + (2*pi);
% Output coordinates
coords = cat(3,flipud(reshape(thOut,size(az))),flipud(reshape(phiOut,size(el))));

% Limit interpolation to points in the selected view
k1 = intersect(find(theta > min(az(:))-0.01),find(theta < max(az(:))+0.01));
k2 = intersect(find(phi > min(el(:))-0.01),find(phi < max(el(:))+0.01));
k = intersect(k1,k2);

% Interpolate to get pixel values at selected points
img = zeros(size(az,1),size(az,2),size(pano,3));
for i = 1:size(pano,3)
    imageValues = pano(:,:,i);
    imageValues = imageValues(:)';
    F = scatteredInterpolant(theta(k)',phi(k)',imageValues(k)','natural');
    %F = TriScatteredInterp(theta(k)',phi(k)',imageValues(k)','linear');
    img(:,:,i) = flipud(reshape(F(az,el),[size(img,1) size(img,2)]));
end

end