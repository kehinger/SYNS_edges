% cleanEdgemapSpheron
%
% Return a list of "good" edge pixels from a Spheron edge map (good = 
% assigned to this view and not outside LiDAR bounds, motion artifact, or
% alignment standard). Note that idx is a set of indices into edges.edge, 
% so edges.edge(idx) are the indices in the image).
%
% 'nocrop' mode keeps pixels not assigned to this view (useful for figures)
% 'nonedge' mode returns locations of "good" pixels (eg, assigned to this
%     view and not artifacts) which are non-edge. If you excluded texture
%     edges before running this code, then the "nonedge" set will include
%     those texture edges.

function idx = cleanEdgemapSpheron(i,edges,maxRep,coords,standards,mode)

lidarBound = -0.25*pi; % lower bound of LiDAR images
lidarMax = 125.5; % values higher than this will be treated as "sky" if not missing

% Identify the pixels that belong to this view (each pixel is assigned
% to the closest camera view; if you skip this step some edges will be
% counted multiple times in overlapping views)
assignedPx = assignPixelsToView(i);
if strcmp(mode,'nocrop')
    assignedPx = ones(size(mask,1),size(mask,2));
end

% Pad the standards map because it's not perfectly accurate
standards = (conv2(double(standards),ones(2),'same') > 0);

% Make separate edge map for each HDR image
edge_map = cell(maxRep,1);
for rep = 1:maxRep
    mag = edges(rep,1).light - edges(rep,1).dark;
    % Find edge pixels at each exposure level
    isEdge = mag > 0;
    % Discard edge pixels not assigned to this view
    badpx = ~assignedPx;
    % Discard edge pixels from the standards
    badpx = (badpx | standards);
    % Discard edge pixels below LiDAR image boundaries
    badpx(coords(:,:,2) < lidarBound) = 1;
    if strcmp(mode,'nonedge')
        edge_map{rep} = find(~isEdge & ~badpx(edges(rep,1).edge));
    else
        edge_map{rep} = find(isEdge & ~badpx(edges(rep,1).edge));
    end
end

% Combine edge maps: use rep1 except where there are artifacts
%if maxRep > 1
%    badpx = single(mask == 1);
%    badpx = (conv2(badpx,ones(artifactDist),'same') > 0);
%    map = badpx.*edge_map(:,:,1) + ((1-badpx).*edge_map(:,:,2));
%else
    map = edge_map{1};
%end

idx = map;

end