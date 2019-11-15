% setExposure
%
% Render an image with a particular exposure and gamma. Output will be
% double in the range 0-1. Input can be any image format.
%
% Input variables:
% im = image
% s = [exposure gamma]
%
% Output variables:
% im = modified image

function im = setExposure(im,s)

im = double(im);
im = (1-exp(-s(1)*im)).^(1/s(2));