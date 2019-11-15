% pcdFillMissing
%
% Returns az, el, and r for each point in the LiDAR image, with the az,el
% coordinates of missing data estimated from the non-missing rows/cols.
%
% Requires the circular statistics toolbox:
% http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
% Input variables:
% pts = point cloud data, size [height width 4]
% varargin = optional flag: 'extend'
%
% Output variables:
% out = matrix of size [height width 4]. Layers are (1) azimuth, (2)
%    elevation, (3) range, and (4) missing data mask (1 if data missing, 0
%    otherwise). Missing data points have range = 200.

function [az,el,r,missing] = pcdFillMissing(pts,varargin)

% Extend elevations
% Some LiDAR images are cropped so their max elevation is less than pi/2
% (since the remaining portion is all sky). Setting the 'extend' flag will
% restore these rows of the LiDAR image with best-guess az and el (and
% misssing r). This is only useful for doing sky/non-sky classification of
% the missing LiDAR values, otherwise this option should be left off.
extendElev = 0;
if ~isempty(varargin)
    if strcmp(varargin{1},'extend')
        extendElev = 1;
    end
end

pts = double(pts);
missing = (pts(:,:,1)==0)&(pts(:,:,2)==0)&(pts(:,:,3)==0)&(pts(:,:,4)==0);

% Convert x,y,z to azimuth and elevation
px = reshape(pts(:,:,1:3),[size(pts,1)*size(pts,2) 3])';
[theta,phi,r] = cart2sph(px(1,:),px(2,:),px(3,:));
theta = reshape(theta,size(pts(:,:,1)));
phi = reshape(phi,size(pts(:,:,1)));
r = reshape(r,size(pts(:,:,1)));

% Use the available points in each row/column to guess the az,el of the
% missing points

% Azimuth
theta_mean = sum(theta)./sum(1-missing);
% Regular averaging will give the wrong answer for columns on the 0/2pi
% seam, so replace these with circular mean (which is slower).
k = find((theta_mean < 0.05*pi)|(theta_mean > 1.95*pi));
for i = k
    theta_mean(1,i) = circ_mean(theta(missing(:,i)==0,i));
end
theta_mean = repmat(theta_mean,[size(pts,1) 1]);

% Elevation
phi_mean = sum(phi,2)./sum(1-missing,2);

% If the 'extend' flag was set and this image is missing higher elevations,
% add rows in the LiDAR image for the missing elevations. The az,el values
% will be filled in with the values estimated from non-missing data. Those
% estimates might be kind of dodgy at the highest elevations because they
% assume no tilt -- it might be better to only fill in a portion of the
% missing elevations (eg, make the output height less than 3771).
if extendElev && length(phi_mean) < 3770
    nRows = 3770 - length(phi_mean);
    phi_mean = vertcat(phi_mean,nan(nRows,1));
    missing = vertcat(missing,ones(nRows,size(missing,2)));
    theta = vertcat(theta,zeros(nRows,size(missing,2)));
    phi = vertcat(phi,zeros(nRows,size(missing,2)));
    r = vertcat(r,zeros(nRows,size(missing,2)));
    theta_mean = vertcat(theta_mean,repmat(theta_mean(end,:),[nRows 1]));
end

% In outdoor scenes, some elevations are entirely missing, but you can
% estimate them by fitting a line to the available elevations.
k = find(~isnan(phi_mean));
if length(k) < length(phi_mean)
    x = 1:length(phi_mean);
    p = polyfit(x(k),phi_mean(k)',1);
    y = (p(1)*x)+p(2);
    n = find(isnan(phi_mean));
    phi_mean(n) = y(n);
end
phi_mean = repmat(phi_mean,[1 size(pts,2)]);

% Reconstructed azimuth, elevation maps
az = (theta.*(1-missing)) + (theta_mean.*missing);
el = (phi.*(1-missing)) + (phi_mean.*missing);
r = (r.*(1-missing)) + (0.*missing);

%{
tic
% Extremely slow when this function is called from another function, but
% fine when this function is used on the command line. Don't know why.
out = zeros(size(pts,1),size(pts,2),4);
out = cat(3,az,el,r,missing);
toc
%}