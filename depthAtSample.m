% depthAtSmaple
%
% Inputs should be size 3 x N

function [d,normD,cornerType,a0] = depthAtSample(p,m1,n1,m2,n2)

% Angle at camera
a0 = acos(dot(normc(m1),normc(m2)));

% Remove edges with no valid surface on one/both sides, and edges where
% the "seed" point was the same on both sides (this can only happen if you
% used a large nearest neighbor search window for fitting the surfaces on
% either side)
%{
badEdge = find((a0==0)|(sum(m1==m2,1)==3)|(sum(m1==0,1)==3)|(sum(m2==0,1)==3));
goodEdge = setdiff(1:size(p,2),badEdge);
%ind = ind(goodEdge);
p = p(:,goodEdge);
m1 = m1(:,goodEdge);
n1 = n1(:,goodEdge);
m2 = m2(:,goodEdge);
n2 = n2(:,goodEdge);
a0 = a0(goodEdge);
%}

% Optional: use the surface means to define the edge plane instead of pixel
% location + theta in the original image
edgeN = cross(m1,m2);
edgeN = edgeN ./ repmat(sqrt(sum(edgeN.^2,1)),[3 1]);

% Line of intersection between edge plane and each surface
v1 = cross(n1,edgeN);
v2 = cross(edgeN,n2);

% Angle at surface
a1 = acos(dot(normc(-m1),v1));
a2 = acos(dot(normc(-m2),v2));

%figure;
%subplot(1,2,1); mask = zeros(imWH); mask(ind) = a1*(180/pi); imagesc(mask); colorbar; axis image off; title('Angle near');
%subplot(1,2,2); mask = zeros(imWH); mask(ind) = a2*(180/pi); imagesc(mask); colorbar; axis image off; title('Angle far');

% Distances to surfaces
dist1 = sqrt(sum(m1.^2));
dist2 = sqrt(sum(m2.^2));

%figure;
%subplot(1,2,1); mask = zeros(imWH); mask(ind) = dist1; imagesc(mask); colorbar; axis image off; title('Dist near');
%subplot(1,2,2); mask = zeros(imWH); mask(ind) = dist2; imagesc(mask); colorbar; axis image off; title('Dist far');

% Expected distances on each side assuming a flat surface extending from
% the other side
b1 = pi-a0-a2; % the angle at the side being estimated
b2 = pi-a0-a1;
b1(b1<=(pi/180)) = pi/180; % require the missing angle to be at least 1 deg (otherwise the missing angle may be tiny or negative due to noisy surface fit)
b2(b2<=(pi/180)) = pi/180;
r1 = dist2.*sin(a2)./sin(b1); % estimated distance on either side
r2 = dist1.*sin(a1)./sin(b2);

%figure;
%subplot(1,2,1); mask = zeros(imWH); mask(ind) = r1; imagesc(mask); colorbar; axis image off; title('Est dist near');
%subplot(1,2,2); mask = zeros(imWH); mask(ind) = r2; imagesc(mask); colorbar; axis image off; title('Est dist far');

% Difference between observed and expected distances
d1 = dist1-r1;
d2 = dist2-r2;

%figure;
%subplot(1,2,1); mask = zeros(imWH); mask(ind) = d1; imagesc(mask); colorbar; axis image off; title('Diff near');
%subplot(1,2,2); mask = zeros(imWH); mask(ind) = d2; imagesc(mask); colorbar; axis image off; title('Diff far');

%disp('Percent of edges detected as corner/curve:');
%sum(sign(d1) == sign(d2))/length(d1)

% If differences have the same sign, the projected surfaces meet in a
% corner or curve. I'll treat these cases as 0 depth change.
% (You also know what kind of corner this is -- if both estimates >
% observed it's a concave corner/curve; if both estimates < observed it's
% a convex corner/curve.)
corner = (sign(d1) == sign(d2))|(d1 == 0)|(d2 == 0);
cornerType = zeros(size(corner));
cornerType((d1<=0)&(d2<=0)) = -1; % concave
cornerType((d1>=0)&(d2>=0)) = 1; % convex

%figure;
%mask = zeros(imWH); mask(ind(concave)) = -1; mask(ind(convex)) = 1; imagesc(mask); axis image off; title('Corners');

%% DEPTH MEASURES

% Basic depth computation: distance between the two surfaces
% This is problematic because it ignores surface slant
d0 = abs(dist1-dist2)';

% Slant-corrected depth: measure depth change at p. This doesn't completely
% remove slant because the impact of measurement noise increases with slant
% It's also problematic because p isn't the true location of the edge, just
% a guess from Spheron.
p1 = dot(m1,n1)./dot(p,n1);
p2 = dot(m2,n2)./dot(p,n2);
%p1 = p - (repmat(dot(p-m1,n1),[3 1]).*n1); % Nope, that's the projection along the normal
%p2 = p - (repmat(dot(p-m2,n2),[3 1]).*n2);
%p1 = sqrt(sum(p1.^2));
%p2 = sqrt(sum(p2.^2));
p_near = p1;
p_far = p2;
p_near(dist1>dist2) = p2(dist1>dist2);
p_far(dist1>dist2) = p1(dist1>dist2);
d = p_far - p_near;

% Slant corrected depth: assume the farther surface is flat and measure
% depth change at the closer surface. This doesn't completely remove slant
% because measurement noise increases with slant, so you are more likely to
% overestimate distance due to noise. (
near_dist = dist1;
near_dist(dist1>dist2) = dist2(dist1>dist2); % closer surface
far_est = r1;
far_est(dist1>dist2) = r2(dist1>dist2); % estimated position of far surface at close surface

%{
% Due to noise, you'll still have issues if both surfaces are very high
% slant (nearly parallel with camera vector).
d(d<0) = 0; % <- remove negative depths -- this probably means flat surface with noise
near_angle = a1;
near_angle(dist1>dist2) = a2(dist1>dist2);
far_angle = b1;
far_angle(dist1>dist2) = b2(dist1>dist2);
normD = repmat(d,[2 1]).*sin(pi-[near_angle; far_angle]); % portion of distance normal to each surface
%d(min(normD) < 0.3) = 0; % smaller normal distance < variance in surface fit -- probably flat surface?
%}

% Change in surface normal
normD = acos(dot(n1,n2))';

%figure;
%subplot(1,2,1); mask = zeros(imWH); mask(ind) = near_angle*(180/pi); imagesc(mask); colorbar; axis image off; title('Angle near');
%subplot(1,2,2); mask = zeros(imWH); mask(ind) = far_angle*(180/pi); imagesc(mask); colorbar; axis image off; title('Angle far');

%figure;
%mask = zeros(imWH,2*imWH); mask(ind) = d0; mask(ind+(1024^2)) = d;
%imagesc(mask); colorbar; axis image off; title('Old depth change vs. slant-corrected depth change');
%colormap(jet);

%figure;
%imagesc(log(mask+1)); colorbar; axis image off; title('Old depth change vs. slant-corrected depth change (log)');
%colormap(jet);


% VARIOUS OPTIONS I DID NOT LIKE
%{
% Old version: compute depth change along the vector projected through the
% edge point. Overestimates depth change on high slant flat surfaces.

% Difference in projection length = separation between planes along sample
% vector. 
p1 = dot((p - m1),n1); % length of projection from p to plane 1
p2 = dot((p - m2),n2); % length of projection from p to plane 2
d = abs(p1-p2);
d(corner) = 0; % depth change is 0 if edges meet in a corner
%d(d>max(d0)) = max(d0); % for figures only

% alternate version because I'm not sure the above is correct...
p1 = dot(m1,n1)./dot(normc(p),n1);
p2 = dot(m2,n2)./dot(normc(p),n2);
d = abs(sqrt(sum(p1.^2,1)) - sqrt(sum(p2.^2,1)));
d(corner) = 0; % depth change is 0 if edges meet in a corner
%d(d>max(d0)) = max(d0); % for figures only

% Nope, that's not right either -- you'd get a scalar, dumbass
% Try this:
% p - (dot(p-m, n) * n)

% Or... average the (observed-expected) distances? This will still
% overestimate depth change on high slant flat surfaces.
d = (abs(d1)+abs(d2))/2;
d(corner) = 0; % depth change is 0 if edges meet in a corner
%d(d>max(d0)) = max(d0); % for figures only
%}

cornerType = cornerType';
a0 = a0';

end


% m1 + k1(v1), m2 + k2(v2) are the two other sides of the quadrilateral
%n1p = n1 - (dot(edgeN,n1)/dot(edgeN,edgeN))*edgeN; % surface normal projected on edge plane
