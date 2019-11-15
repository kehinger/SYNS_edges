% edgeCharacterization
%
% Characterize edges in Spheron or LiDAR images by fitting surfaces on
% either side of edge.
%
% To run a list of edges, pass the view number and indices in the image:
%     edgeCharacterization(SYNSno,viewNo,edgeIndices)
%
% Otherwise, this function will process all edges from all views of a
% specified scene:
%     edgeCharacterization(SYNSno)

function [surfaceData,p,seedPts] = edgeCharacterization(folderNo,varargin)

windowPx = [32 32]; % base size of LiDAR patch used in surface fit
seedOffset = 0.0024; % move this distance (in rad) from edge before starting surface fit on each side
closestSampleThresh = 0.0009; % max distance for "seed" point on each side of edge
outlierThresh = 0.5; % proportion of outliers allowed in the convex hull of a "good" surface
maxK = 300; % neighborhood size for surface fit
skyDist = 9999;
folderName = 'ground_truth_nz3_step2';
imageSource = 'Spheron';
edgesToRun = [];

load excludedViews

% To run on a single view and/or list of matched LiDAR edges:
if ~isempty(varargin)
    viewsToRun = varargin{1};
    modeSpheron = 'nocrop';
    if length(varargin) > 1
        edgesToRun = varargin{2};
    	imageSource = 'LiDAR';
    end
% DEFAULT: run on in Spheron mode on all available views
else
    modeSpheron = 'standard';
    viewsToRun = setdiff(1:42,[27 excludedViews{folderNo}]);
end

% Load LiDAR point cloud
pts = loadpcdFixed(['SYNS/' int2str(folderNo) '/clean_cloud.pcd']);
nonMissing = find(sum(abs(pts(:,1,:)),3)>0);
% Apply manual co-regsitration correction
if exist(['SYNS/' int2str(folderNo) '/shift_proposed.mat'],'file')
    disp('Applying manual co-registration correction');
    proposed_shift = load(['SYNS/' int2str(folderNo) '/shift_proposed.mat']);
    for ch = 1:3
        coords = pts(:,:,ch);
        coords(nonMissing) = coords(nonMissing)+proposed_shift.s1(ch);
        pts(:,:,ch) = coords;
    end
    clear proposed_shift
else
    disp('No manual co-registration correction');
end

[offset,~,~] = cart2sph(pts(nonMissing,1,1),pts(nonMissing,1,2),pts(nonMissing,1,3));
offset = mean(offset); % azimuth of left edge of LiDAR image
originalIndices = 1:(size(pts,1)*size(pts,2));
originalIndices = reshape(originalIndices(:),[size(pts,1) size(pts,2)]);

% Fill in sky pixels (outdoor images only)
if folderNo <= 80
    % Load sky mask
    skymask = imread(['SYNS/' int2str(folderNo) '/sky_mask.png']);
    skymask = skymask(1:size(pts,1),:);
    skymask = (skymask > 200);
    
    % Fill in approximate locations for the sky pixels and mark them as
    % extremely high depth
    [az,el,~,missing] = pcdFillMissing(pts);
    issky = (missing & (skymask == 1));
    [x,y,z] = sph2cart(az(issky),el(issky),skyDist*(ones(sum(issky(:)),1)));
    skyXYZ = [x y z];
    for i = 1:3
        panel = pts(:,:,i);
        panel(issky) = skyXYZ(:,i);
        pts(:,:,i) = panel;
    end
end

if ~exist(['SYNS/' int2str(folderNo) '/' folderName],'dir')
    mkdir(['SYNS/' int2str(folderNo) '/' folderName]);
end

for viewNo = viewsToRun
if ~exist(['SYNS/' int2str(folderNo) '/' folderName '/edgesElder' int2str(viewNo) '.mat'],'file')

    % LiDAR version
    if strcmp(imageSource,'LiDAR')
        % Load view and edge data
        load(['projectionLiDAR/view' int2str(viewNo) '.mat'],'coords','view');
        load(['SYNS/' int2str(folderNo) '/lidar_projections/view' int2str(viewNo) '.mat'],'img','alpha','sky');
        load(['SYNS/' int2str(folderNo) '/lidar_projections/edges' int2str(viewNo) '.mat'],'edges');
        load(['SYNS/' int2str(folderNo) '/projected_views/view' int2str(viewNo) '.mat'],'standards');
        imWH = size(img,1);
        
        % Clean up the set of LiDAR edges
        if isempty(edgesToRun)
            [idx,~] = cleanEdgemapLiDAR(viewNo,edges,coords,img,alpha,standards,sky,'standard');
        else
            [~,idx] = ismember(edgesToRun,edges.edge);
        end
        clear img alpha sky standards coords
        
    % Spheron version
    elseif strcmp(imageSource,'Spheron')
        % Load view and edge data
        load(['projectionSpheron/view' int2str(viewNo) '.mat'],'coords','view');
        load(['SYNS/' int2str(folderNo) '/projected_views/view' int2str(viewNo) '.mat'],'img','standards');
        load(['SYNS/' int2str(folderNo) '/projected_views/edges_ElderZucker3_' int2str(viewNo) '.mat']);
        imWH = size(img,1);
        
        % Clean Spheron edge map
        edges = edges(1,1);
        idx = cleanEdgemapSpheron(viewNo,edges,1,coords,standards,modeSpheron);
        clear img standards coords
        
    end
    
    % Edge info
    lightdark = [edges.light(idx) edges.dark(idx)]; % lightdark depth
    theta = edges.g1dir(idx); % angle
    subpx = [edges.xzero(idx) edges.yzero(idx)]; % subpixel (x,y) edge locations
    
    %% SURFACE FIT
    % Edge points in LiDAR xyz coordinates
    p = subpixel2xyz(subpx-0.5,[imWH imWH],view)';
    p(2,:) = -p(2,:);
    
    % Adjust LiDAR window size for elevation (larger near poles)
    [az,el,~] = cart2sph(p(1,:),p(2,:),p(3,:));
    adj = round(windowPx(2)*sin(el/2).^2);
    
    % Get LiDAR image x,y for each edge point's az,el
    lidary = round((el*(1258/.7864)) + 1259);
    lidarx = mod(round((offset-az)*(5027/pi)),10054)+1;
    
    % Edge normals (= normal of a plane that includes all points along the edge
    % plus the camera at (0,0,0))
    q = subpixel2xyz(subpx-0.5+[sin(theta) cos(theta)],[imWH imWH],view)';
    q(2,:) = -q(2,:);
    edgeN = cross(p,q)';
    edgeN = edgeN ./ repmat(sqrt(sum(edgeN.^2,2)),[1 3]);
    
    % Run surface fit at each edge point
    surfaceData = zeros(length(idx),14);
    seedPts = zeros(length(idx),10);
    inliers = zeros(length(idx),8);
    surfaceDataExtended = zeros(length(idx),8);
    surfaceChecks = zeros(length(idx),4);
    for j = 1:length(idx)
    %while 1
        %j = randi(length(idx),1);
        
        % Get LiDAR data in a rectangular patch around edge pixel.
        % Searching for nearest neighbor in this patch is much, much faster
        % than searching the entire LiDAR image.
        win = windowPx + [0 adj(j)]; % adjusted window
        if (lidary(j)-win(1)) < 1
            a = 1:2*win(1);
        elseif (lidary(j)+win(1)) > size(pts,1)
            a = (size(pts,1)-(2*win(1))):size(pts,1);
        else
            a = (lidary(j)-win(1)):(lidary(j)+(win(1)));
        end
        if (lidarx(j)-win(2)) < 1
            b = [1:(lidarx(j)+win(2)) (size(pts,2)+lidarx(j)-win(2)):size(pts,2)];
        elseif (lidarx(j)+win(2)) > size(pts,2)
            b = [1:(lidarx(j)+win(2)-size(pts,2)) (lidarx(j)-win(2)):size(pts,2)];
        else
            b = (lidarx(j)-win(2)):(lidarx(j)+(win(2)));
        end
        patch_orig = pts(a,b,4);
        patch_ind_orig = 1:(size(patch_orig,1)*size(patch_orig,2));
        patch = reshape(pts(a,b,1:3),[length(a)*length(b) 3]);
        patchIndices = reshape(originalIndices(a,b),[length(a)*length(b) 1]);
        
        % Drop missing points
        keepPts = (sum(abs(patch),2) > 0);
        patch = double(patch(keepPts,:));
        patchIndices = patchIndices(keepPts,:);
        patch_ind_orig = patch_ind_orig(keepPts);
        
        %edgeLabel = sign(dot(patch'-repmat(p(:,j),[1 size(patch,1)]),repmat(edgeN(j,:)',[1 size(patch,1)])));

        % Fit surfaces on either side of edge
        if size(patch,1) >= maxK
            
            % Seed the surface fit with points offset from the edge
            % Use whichever LiDAR point is closest to the desired seed position
            p1 = p(:,j) + (seedOffset*edgeN(j,:)');
            p2 = p(:,j) - (seedOffset*edgeN(j,:)');
            d1 = abs(acos(dot(normc(patch'),repmat(p1/norm(p1),[1 size(patch,1)]))));
            d2 = abs(acos(dot(normc(patch'),repmat(p2/norm(p2),[1 size(patch,1)]))));
            [~,ord1] = sort(d1);
            [~,ord2] = sort(d2);
            seed1 = patch(ord1(1),:); % position of closest valid LiDAR point
            seed2 = patch(ord2(1),:);
            seedDist1 = d1(ord1(1)); % angular distance to closest valid LiDAR point
            seedDist2 = d2(ord2(1));
            %seed1_for_plot = patch_ind_orig(ord1(1));
            %seed2_for_plot = patch_ind_orig(ord2(1));
            seedPts(j,:) = [patchIndices(ord1(1)) seedDist1 seed1 patchIndices(ord2(1)) seedDist2 seed2];
            
            % Order patch points by distance to each seed and run surface fit
            % Use spherical distance and DO NOT restrict the set to points on
            % that side of edge
            if (seedDist1 < closestSampleThresh)
                d1 = sum((repmat(seed1,[size(patch,1) 1])-patch).^2,2);
                [~,ord] = sort(d1,'ascend');
                [knn,fit_norm,fit_mean,inK] = planefit(double(seed1(1)),double(seed1(2)),double(seed1(3)),...
                    double(patch(ord(1:300),1)),double(patch(ord(1:300),2)),double(patch(ord(1:300),3)));
                if sum((fit_mean+fit_norm).^2) > sum((fit_mean-fit_norm).^2)
                    fit_norm = -fit_norm;
                end
                plane1 = [knn fit_mean fit_norm];
                inK1 = ord(inK==1);
            else
                plane1 = zeros(1,7);
                inK1 = [];
            end
            if (seedDist2 < closestSampleThresh)
                d2 = sum((repmat(seed2,[size(patch,1) 1])-patch).^2,2);
                [~,ord] = sort(d2,'ascend');
                [knn,fit_norm,fit_mean,inK] = planefit(double(seed2(1)),double(seed2(2)),double(seed2(3)),...
                    double(patch(ord(1:300),1)),double(patch(ord(1:300),2)),double(patch(ord(1:300),3)));
                if sum((fit_mean+fit_norm).^2) > sum((fit_mean-fit_norm).^2)
                    fit_norm = -fit_norm;
                end
                plane2 = [knn fit_mean fit_norm];
                inK2 = ord(inK==1);
            else
                plane2 = zeros(1,7);
                inK2 = [];
            end
            
            % Reject surfaces where all points are colinear in the LiDAR
            % map (same row or same column of pixels). There was a check
            % for this in Alex's code, but I omitted it because my version
            % uses non-grid data and has no way of knowing if some subset
            % of samples were from the same pixel row/column
            [pts_y, pts_x] = ind2sub(size(patch_orig),patch_ind_orig);
            if (length(unique(pts_y(inK1))) == 1)||(length(unique(pts_x(inK1))) == 1)
                inK1 = [];
                plane1 = zeros(1,7);
            end
            if (length(unique(pts_y(inK2))) == 1)||(length(unique(pts_x(inK2))) == 1)
                inK2 = [];
                plane2 = zeros(1,7);
            end
            
            % Reject surfaces which are very sparse or oddly-shaped --
            % these occur where there is a complex occluder like chicken
            % wire, or where the code collects co-planar points from many
            % different objects. Depth estimates from these surfaces are
            % very unreliable.
            if ~isempty(inK1)
                %pts_xy = points2plane(patch',m1(:,1),n1(:,1))'; % project patch onto surface 1
                chull = convhull(pts_y(inK1),pts_x(inK1)); % convex hull of surface
                outliers = inpolygon(pts_y(setdiff(1:size(patch,1),inK1)),pts_x(setdiff(1:size(patch,1),inK1)),...
                    pts_y(inK1(chull)),pts_x(inK1(chull))); % non-member points in the convex hull
                if (sum(outliers)/length(inK1)) > outlierThresh
                    inK1 = [];
                    plane1(1) = -1;
                end
                surfaceChecks(j,1) = sum(outliers);
            end
            if ~isempty(inK2)
                %pts_xy = points2plane(patch',m2(:,1),n2(:,1))'; % project patch onto surface 2
                chull = convhull(pts_y(inK2),pts_x(inK2)); % convex hull of surface
                outliers = inpolygon(pts_y(setdiff(1:size(patch,1),inK2)),pts_x(setdiff(1:size(patch,1),inK2)),...
                    pts_y(inK2(chull)),pts_x(inK2(chull))); % non-member points in the convex hull
                if (sum(outliers)/length(inK2)) > outlierThresh
                    inK2 = [];
                    plane1(1) = -1;
                end
                surfaceChecks(j,2) = sum(outliers);
            end
            
            surfaceData(j,:) = [plane1 plane2];
            
            % Identify the points on either side of edge
            plane1_pts = patch(inK1,:);
            %plane1_pts_ind = patch_ind_orig(inK1);
            plane2_pts = patch(inK2,:);
            %plane2_pts_ind = patch_ind_orig(inK2);
            
            % If surface fit succeeds on both sides, do further steps to
            % characterize the neighborhood around the edge
            if (sum(abs(plane1(2:end)))>0) && (sum(abs(plane2(2:end)))>0) && (~isempty(inK1)) && (~isempty(inK2))
                
                % Compute standard deviation of inliers on each plane
                prj_dist1 = dot(repmat(plane1(5:7),[size(plane1_pts,1) 1])',(plane1_pts-repmat(plane1(2:4),[size(plane1_pts,1) 1]))');
                std1 = sqrt(sum(prj_dist1.^2)/(length(prj_dist1) - 3));
                prj_dist2 = dot(repmat(plane2(5:7),[size(plane2_pts,1) 1])',(plane2_pts-repmat(plane2(2:4),[size(plane2_pts,1) 1]))');
                std2 = sqrt(sum(prj_dist2.^2)/(length(prj_dist2) - 3));
                
                surfaceChecks(j,3:4) = [std1 std2];
                
                % Identify a set of points in between the two surfaces. I
                % do this by placing a disk between the two means (diameter
                % = distance between the means).
                ctr = (normc(plane1(2:4)')+normc(plane2(2:4)'))/2;
                angTh = acos(dot(ctr,normc(plane1(2:4)')));
                ang = acos(dot((patch./repmat(sqrt(sum(patch.^2,2)),[1 3]))',repmat(ctr,[1 size(patch,1)])));
                extra = setdiff(find(ang<angTh),union(inK1,inK2));
                extra_pts = patch(extra,:);
                
                % Assign the in-between points to the two planes
                prj_dist1_extended = dot(repmat(plane1(5:7),[size(extra_pts,1) 1])',(extra_pts-repmat(plane1(2:4),[size(extra_pts,1) 1]))');
                prj_dist2_extended = dot(repmat(plane2(5:7),[size(extra_pts,1) 1])',(extra_pts-repmat(plane2(2:4),[size(extra_pts,1) 1]))');
                inK1_extended = union(inK1,extra(abs(prj_dist1_extended) < 3*std1));
                inK2_extended = union(inK2,extra(abs(prj_dist2_extended) < 3*std2));
                %unassigned = setdiff(extra,union(inK1_extended,inK2_extended));
                plane1_pts_extended = patch(inK1_extended,:);
                %plane1_pts_extended_ind = patch_ind_orig(inK1_extended);
                plane2_pts_extended = patch(inK2_extended,:);
                %plane2_pts_extended_ind = patch_ind_orig(inK2_extended);
                
                % Save data about the extended surface
                surfaceDataExtended(j,:) = [length(inK1_extended) mean(plane1_pts_extended) length(inK2_extended) mean(plane2_pts_extended)];
                
                % Inlier check: check how many of surface 1's K inlier
                % points would be inliers on surface 2 and vice versa
                % Compute projection distances for pts1 onto plane2 and vice versa
                dist_1on2 = dot(repmat(plane2(5:7),[size(plane1_pts,1) 1])',(plane1_pts-repmat(plane2(2:4),[size(plane1_pts,1) 1]))');
                dist_2on1 = dot(repmat(plane1(5:7),[size(plane2_pts,1) 1])',(plane2_pts-repmat(plane1(2:4),[size(plane2_pts,1) 1]))');
                % How many points would be inliers on the other plane?
                inliers(j,1) = sum(abs(dist_1on2) < 3*std2);
                inliers(j,2) = sum(abs(dist_2on1) < 3*std1);
                % T-test: projection to own plane vs. projection to other plane
                [~,pval,~,~] = ttest2(prj_dist1',dist_2on1');
                inliers(j,3) = pval;
                [~,pval,~,~] = ttest2(prj_dist2',dist_1on2');
                inliers(j,4) = pval;
                
                % Inlier check using extended neighborhood
                % Compute projection distances for pts1 onto plane2 and vice versa
                dist_1on2_extended = dot(repmat(plane2(5:7),[size(plane1_pts_extended,1) 1])',(plane1_pts_extended-repmat(plane2(2:4),[size(plane1_pts_extended,1) 1]))');
                dist_2on1_extended = dot(repmat(plane1(5:7),[size(plane2_pts_extended,1) 1])',(plane2_pts_extended-repmat(plane1(2:4),[size(plane2_pts_extended,1) 1]))');
                % How many points would be inliers on the other plane?
                inliers(j,5) = sum(abs(dist_1on2_extended) < 3*std2);
                inliers(j,6) = sum(abs(dist_2on1_extended) < 3*std1);
                % T-test: projection to own plane vs. projection to other plane
                [~,pval,~,~] = ttest2(prj_dist1_extended',dist_2on1_extended');
                inliers(j,7) = pval;
                [~,pval,~,~] = ttest2(prj_dist2_extended',dist_1on2_extended');
                inliers(j,8) = pval;
                
            end
        end
    end
    p = p';
    edgeN = edgeN';
    if isempty(edgesToRun)
        save(['SYNS/' int2str(folderNo) '/' folderName '/edgesElder' int2str(viewNo) '.mat'],'idx','lightdark','theta','subpx','p','edgeN','seedPts','surfaceData','surfaceChecks','surfaceDataExtended','inliers');
    end
end
end

end

% Convert subpixel locations in image to LiDAR coordinates
% pxy = N x 2 list of (x,y) values
% imSize = [imageWidth imageHeight]
function p = subpixel2xyz(pxy,imSize,view)
	view.azel = (pi/180)*view.azel;
    %view.azel(1) = -view.azel(1); % reverse az
    halfFOV = (pi/180)*(view.fov/2);
    pxy = (pxy-repmat(imSize/2,[size(pxy,1) 1])).*repmat((2./imSize),[size(pxy,1) 1]);
    pxy(:,2) = -pxy(:,2); % flip up/down because +el means up
    pxy = pxy.*repmat(tan(halfFOV),[size(pxy,1) 1]);
    p = rotElAz([ones(size(pxy,1),1) pxy(:,1) pxy(:,2)]',-view.azel)';
    p = p./repmat(sqrt(sum(p.^2,2)),[1 size(p,2)]);
end

% Return x,y positions of points projected to a 2D plane. The mean of the
% plane (m) is used as the centre (0,0) of this plane and the positive x
% axis points towards 0 in the original, 3D coordinate frame.
function xy = points2plane(pts,m,n)
    p1 = pts-(dot((pts-repmat(m,[1 size(pts,2)])),repmat(n,[1 size(pts,2)])).*repmat(n,[1 size(pts,2)]));
    p2 = -m - (dot(-m,n)*n);
    r = sqrt(sum((p1-repmat(p2,[1 size(p1,2)])).^2,1));
    th = acos(dot(normc(p1),repmat(normc(p2),[1 size(p1,2)])));
    thcross = cross(normc(p1),repmat(normc(p2),[1 size(p1,2)]));
    thsign = sign(dot(repmat(n,[1 size(pts,2)]),thcross));
    th = th.*thsign;
    xy = [r.*cos(th); r.*sin(th)];
end
