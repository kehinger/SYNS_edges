function demo

% This demo runs the entire processing pipeline on a single SYNS scene
% (takes several hours).

SYNSno = 1; % which scene to run

% Project a set of standard camera views from the Spheron image and detect
% edges in these views.
detectEdges(SYNSno);

% Locate the Spheron edge in the LiDAR point cloud and try to fit a plane
% on either side of the edge.
edgeCharacterization(SYNSno);

% Using the plane fits, label edges as depth vs. non-depth. Exclude bad
% edges and compute depth contrast.
compute_depth_contrast_with_exclusions(SYNSno);

end
