% pixelCoords
%
% Gets azimuth,elevation coordinates for pixels in a camera view through a
% spherical panorama. The resolution of the camera view is automatically
% determined by the number of pixels in the panorama.
%
% Input variables:
% panosize = [H W C], size of the panoramic image in equirectangular
%   projection
% fov = [horz vert] specifies FOV of the camera (in degrees)
%
% Output variables:
% az,el = coordinates of pixels in the projected camera view

function [az,el] = pixelCoords(panosize,fov,outSize)

% Number of sample points in view (= pixels in output image)
% Note in case you want to compute output size automatically:
% The limiting factor is the number of pixels in the left/right column of a
% view at the equator (this is always smaller than the number of pixels in
% the center of the view, and goes to 0 as horizontal f.o.v. goes to 180).
fW = (pi/180)*(fov(1)/2);
fH = (pi/180)*(fov(2)/2);
minVertPx = outSize(1); % use the user-specified output size
imgH = round(2*minVertPx);
imgW = round(imgH*(fov(1)/fov(2)));

% Segments of viewing plane (correspond to pixels in output image):
s1 = ((1:imgW)-0.5-(imgW/2))*(2/imgW);
s2 = ((1:imgH)-0.5-(imgH/2))*(2/imgH);

% y,z positions of segements (image plane is x=1)
s1 = s1*tan(fW);
s2 = s2*tan(fH);
s1 = repmat(s1,[imgH 1]);
s2 = repmat(s2',[1 imgW]);

[az,el,~] = cart2sph(ones(imgH*imgW,1),s1(:),s2(:));
az = reshape(az,[imgH imgW]);
el = reshape(el,[imgH imgW]);

end