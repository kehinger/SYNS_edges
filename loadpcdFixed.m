% loadpcdFixed
%
% Loads a point cloud (.pcd file) and reshapes the output to fix an issue
% with the default version of loadpcd.
%
% Requires matpcl functions available here:
% https://www.mathworks.com/matlabcentral/fileexchange/40382-matlab-to-point-cloud-library

function pts = loadpcdFixed(filename)

pts = loadpcd(filename);
% The point cloud data is the right size (~3000 x ~10000) but the data is
% read in in the wrong orientation -- columns are slotted into the matrix
% as rows. Not sure if this is caused by the .pcd file headers/format or a
% bug in loadpcd.

[h,w,c] = size(pts);
pts = permute(pts,[2 1 3]); % transpose height/width
pts = reshape(pts,[h*w c]); % turn into columns
pts = reshape(pts,[h w c]); % and back to height x width

% Remove any data in the first row of the LiDAR image (any non-missing
% values in this row are noise)
pts(1,:,:) = 0;

end
