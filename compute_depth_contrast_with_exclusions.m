% compute_depth_contrast_with_exclusions
%
% Post-process the edge data files produced by edgeCharacterization. Detect
% edge type (non-depth vs depth) and compute depth difference across edge.
% Exclude bad edges.

function compute_depth_contrast_with_exclusions(i)

load excludedViews
skyTh = 200; % surfaces beyond this distance are assumed to be sky
hasArtifacts = 0;
if exist(['SYNS/' int2str(i) '/projected_views/rep1_artifacts.mat'],'file')
    hasArtifacts = 1;
    artifacts = load(['SYNS/' int2str(i) '/projected_views/rep1_artifacts.mat']);
end

%% Combine the edge data files for all views in this scene
% Also compute depth/normal change and do some quality checks to exclude
% bad edges.
for j = setdiff(1:42,[27 excludedViews{i}])
    
    load(['SYNS/' int2str(i) '/ground_truth_nz3_step2/edgesElder' int2str(j) '.mat']);
    load(['SYNS/' int2str(i) '/projected_views/edges_ElderZucker3_' int2str(j) '.mat']);
    
    rep1edges = edges(1,1);
    scene = i*ones(size(idx));
    view = j*ones(size(idx));
    edge = rep1edges.edge(idx);
    rep = ones(size(idx));
    [y,x] = ind2sub([640 640],edge);
    xy = [x y];
    coords = p;
    [az,el,~] = cart2sph(coords(:,1),coords(:,2),coords(:,3));
    azel = [az el];
    edgeN = edgeN';
    seedInd = seedPts(:,[1 6]);
    seedDist = seedPts(:,[2 7]);
    seed1 = seedPts(:,3:5);
    seed2 = seedPts(:,8:10);
    knn = surfaceData(:,[1 8]);
    m1 = surfaceData(:,2:4);
    m2 = surfaceData(:,9:11);
    n1 = surfaceData(:,5:7);
    n2 = surfaceData(:,12:14);
    r = [sqrt(sum(m1.^2,2)) sqrt(sum(m2.^2,2))];
    
    % Find cases where surfaces fall on the same side of Spheron edge
    sideOfEdge1 = sign(dot(m1'-coords',edgeN'))';
    sideOfEdge2 = sign(dot(m2'-coords',edgeN'))';
    sameSide = (sideOfEdge1 == sideOfEdge2);
    
    % Extended surface data
    knnExt = surfaceDataExtended(:,[1 5]);
    m1Ext = surfaceDataExtended(:,2:4);
    m2Ext = surfaceDataExtended(:,6:8);
    rExt = [sqrt(sum(m1Ext.^2,2)) sqrt(sum(m2Ext.^2,2))];
    
    % Variables used to exclude bad edges
    chullOutliers = surfaceChecks(:,1:2);
    projStdDev = surfaceChecks(:,3:4);
    ttestPval = inliers(:,[3:4 7:8]);
    oneSide = (r(:,1) <= 0)|(r(:,2) <= 0)|(knn(:,1) < 1)|(knn(:,2) < 1)|sameSide; % oneSide if either surface is missing or on wrong side
    skyEdge = (r(:,1) > skyTh)|(r(:,2) > skyTh); % sky if either distance > some threshold
    
    % Extra exclusions
    seedDistExclusion = (min(seedDist,[],2) >= .0009); % reject edge if no LiDAR samples within 1.5x width of a LiDAR pixel on one side
    bothSky = (r(:,1) > skyTh)&(r(:,2) > skyTh); % reject edges if both sides are sky
    oneSide = oneSide | seedDistExclusion | bothSky;
    
    % Inliers check for depth edge
    % There is no advantage to using the extended KNN for this check, so I
    % use the original KNN
    inliers = inliers(:,[1:2 5:6]);
    inliersPerc = inliers./[knn knnExt];
    isDepth = (sum(inliers(:,1:2),2) == 0); % original neighborhoods
    
    % Depth and normal change variables
    % There is no advantage to using the extended KNN for the corner check,
    % so I use the original KNN
    depthChange = abs(r(:,1)-r(:,2));
    [~,normalChange,corner,camAngle] = depthAtSample(p',m1',n1',m2',n2');
    depthChange(corner ~= 0) = 0;
    depthChange(~isDepth) = 0; % depth = 0 for edges that share inliers OR are corners
    
    % Depth change from extended KNN is probably slightly more accurate --
    % use this instead of the other "depthChange" variable
    depthChangeExt = abs(rExt(:,1)-rExt(:,2));
    depthChangeExt(corner ~= 0) = 0;
    depthChangeExt(~isDepth) = 0; % depth = 0 for edges that share inliers OR are corners
    
    newT = table(scene,view,idx,edge,rep,xy,subpx,theta,lightdark,coords,azel,oneSide,skyEdge,edgeN,seedInd,seedDist,seed1,seed2,chullOutliers,projStdDev,ttestPval,inliers,inliersPerc,knn,m1,m2,n1,n2,r,isDepth,corner,depthChange,normalChange,camAngle,knnExt,m1Ext,m2Ext,rExt,depthChangeExt);
    
    % Remove edges that fall within the artifact mask (if any)
    if hasArtifacts
        toRemove = artifacts.maskInd{j};
        if ~isempty(toRemove)
            goodInd = ~ismember(newT.edge,toRemove);
            newT = newT(goodInd,:);
        end
    end
    
    if exist('T','var')
        T = vertcat(T,newT);
    else
        T = newT;
    end
    
end

% Save the cleaned-up ground truth file
save(['SYNS/' int2str(i) '/ground_truth_nz3_step2.mat'],'T');
s = size(T,1); % Number of edges before "foliage" exclusions

%% Compute ground truth measure (depth contrast) and exclude foliage edges

% Exclude bad edge fits (mostly foliage). Edges are excluded if there are
% too few inliers on the smaller side, too far away, or have too many
% non-inlier points within either surface's convex hull.
chullTest = max(T.chullOutliers./(T.chullOutliers+T.knn),[],2);
toKeep = (min(T.rExt,[],2) <= 50) & (chullTest < .15) & (min(T.knn,[],2) > 14);
T = T(toKeep,:);

% Compute depth contrast using extended KNN, set corners and edges with
% common inlier points to "non-depth"
T.depthContrast = T.depthChangeExt ./ sum(T.rExt,2);
T.depthContrast(T.isDepth == 0) = 0;
T.depthContrast(T.corner ~= 0) = 0;

% Add the figure-ground label
% NOTE: seed1 is on the dark side of the detected edge and seed2 is on the
% light side. lightNear = 1 means light side is closer; lightNear = -1
% means dark side is closer.
T.lightNear = sign(T.rExt(:,1)-T.rExt(:,2)); % is the dark side farther / light side closer?

% To speed up later processing, remove failed surface fits
T = T(T.oneSide == 0,:);

save(['SYNS/' int2str(i) '/ground_truth_nz3_exclusions.mat'],'T');
f = size(T,1); % Number of edges after "foliage" exclusions
disp(['Scene ' int2str(i) ': ' int2str(s) ' -> ' int2str(f) ' edges after exclusions']);

end
