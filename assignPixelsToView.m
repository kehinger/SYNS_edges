function mask = assignPixelsToView(i)
% Make a mask showing which pixels below to this view.

s = sphere_icos2_points(2); % 42 views
c = load(['projectionSpheron/view' int2str(i) '.mat']);
az = c.coords(:,:,1); az = az(:);
el = c.coords(:,:,2); el = el(:);
[x,y,z] = sph2cart(az,el,ones(size(az)));
px = normc([x y z]');
angDist = zeros(size(s,2),size(px,2));
for j = 1:size(s,2)
    d = abs(acos(dot(px,repmat(s(:,j),[1 size(px,2)]))));
    angDist(j,:) = d;
end
[~,closest] = min(angDist);
closest = reshape(closest,[size(c.coords,1) size(c.coords,2)]);
mask = (closest == i);