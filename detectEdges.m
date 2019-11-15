% detectEdges
%
% Detect edges in projected views of a SYNS Spheron image. If the projected
% views does not exist for this scene, this function will create them.

function detectEdges(folderNo)

% List of views to skip (no valid LiDAR returns)
load excludedViews

%% 1. Project ~40 standard camera views from each Spheron image
% Views are 48 x 48 degrees, projected to 640 x 640 px images.
cameraFOV = 48; % degrees
viewWH = 640;
if isempty(dir(['SYNS/' int2str(folderNo) '/projected_views/view*.mat']))
    
    if ~exist(['SYNS/' int2str(folderNo) '/projected_views'],'dir')
        mkdir(['SYNS/' int2str(folderNo) '/projected_views']);
    end
    fov = [cameraFOV cameraFOV];
    
    % Sample evenly over sphere
    s = sphere_icos2_points(2); % gives 42 vectors w/ ~30 degrees between neighboring views
    
    % Convert x,y,z to az,el (for pano2view function)
    [az,el,~] = cart2sph(s(1,:)',s(2,:)',s(3,:)');
    el = (180/pi)*el; az = (180/pi)*az;
    
    % Read the HDR images (note that some scenes do not have a rep2.hdr)
    a = hdrread(['SYNS/' int2str(folderNo) '/rep1.hdr']);
    try
        b = hdrread(['SYNS/' int2str(folderNo) '/rep2.hdr']);
        rep2exists = 1;
    catch
        rep2exists = 0;
    end
    
    % Load the hand-annotated masks. Every image has a standards mask which
    % removes the 3 calibration targets+tripods. Some images also have a
    % mask for the rep1 HDR image which removes movement artifacts and
    % objects which moved in between the LiDAR and Spheron image capture.
    maskStandards = double(imread(['SYNS/' int2str(folderNo) '/standards_mask.png']));
    maskStandards = (conv2(maskStandards,ones(2),'same') > 0); % extend the mask a little to remove background pixels near standards
    try
        maskImg = double(imread(['SYNS/' int2str(folderNo) '/rep1_mask.png']));
        maskImg = (maskImg > 0);
        maskexists = 1;
    catch
        maskexists = 0;
    end
    
    % Project views, save HDR images, projected masks, and az,el
    % coordinates
    for i = setdiff(1:42,[27 excludedViews{folderNo}])
        
        % Project views (rep1 and 2 if available)
        [img1,coords] = pano2view(double(a),[az(i) el(i)],fov,viewWH);
        if rep2exists
            [img2,~] = pano2view(double(b),[az(i) el(i)],fov,viewWH);
            img = cat(3,img1,img2);
        else
            img = img1;
        end
        view.vector = s(:,i);
        view.azel = [az(i) el(i)];
        view.fov = fov;
        
        % Project standards mask
        [standards,~] = pano2view(double(maskStandards),[az(i) el(i)],fov,viewWH);
        standards = (standards > 0.5);
        
        % Project artifacts mask, if it exists, or create an empty mask
        if maskexists
            [artifacts,~] = pano2view(double(maskImg),[az(i) el(i)],fov,viewWH);
            artifacts = (artifacts > 0.5);
            artifacts = (conv2(double(artifacts),ones(2),'same') > 0); % add some padding to exclude more pixels near artifacts
        else
            artifacts = false(viewWH);
        end
        
        % Check if any of the pixels assigned to this view contain
        % artifacts, and save their indices
        assignedPx = assignPixelsToView(i);
        artifactsInView = (artifacts & assignedPx);
        if sum(artifactsInView(:)) > 0
            maskInd{i} = find(artifacts == 1);
        else
            maskInd{i} = [];
        end
        
        % Save the image(s) and view details
        save(['SYNS/' int2str(folderNo) '/projected_views/view' int2str(i) '.mat'],'img','standards','view');
        imwrite(setExposure(img1,[5 1]),['SYNS/' int2str(folderNo) '/projected_views/im' int2str(i) '-1.png']);
        if rep2exists
            imwrite(setExposure(img2,[5 1]),['SYNS/' int2str(folderNo) '/projected_views/im' int2str(i) '-2.png']);
        end
        
        % Save the coordinates of individual pixels in the view for later
        % reference. Note that the coordinates for a given view number
        % (e.g., 16) are identical for all spheres, so you only need to
        % save them once.
        if ~exist(['projectionSpheron/view' int2str(i) '.mat'],'file')
            if ~exist('projectionSpheron','dir')
                mkdir('projectionSpheron');
            end
            save(['projectionSpheron/view' int2str(i) '.mat'],'coords','view');
        end
        
        clear img1 img coords view standards
        if rep2exists
            clear img2
        end
        
    end
    
    % Save the artifact masks for all views in a single file. These masks
    % will be used to exclude edges after combining edge detection results
    % from all views, so it's more convenient to have them in a single
    % file.
    save(['SYNS/' int2str(folderNo) '/projected_views/rep1_artifacts.mat'],'maskInd');

end

%% 2. Detect edges in each projected view
% Edges are detected with the Elder-Zucker algorithm with noise=3.
addpath(genpath('edgeElder'));
nz = 3;
if isempty(dir(['SYNS/' int2str(folderNo) '/projected_views/edges_ElderZucker' int2str(nz) '_*.mat']))
    
    % Detect edges in all views
    for i = setdiff(1:42,[27 excludedViews{folderNo}])
        
        % Load projected view
        load(['SYNS/' int2str(folderNo) '/projected_views/view' int2str(i) '.mat']);
        maxRep = 1;
        if size(img,3) == 6
            maxRep = 2;
        end
        
        % Run edge detector on HDR image (convert to uint8 so that median
        % pixel value = 0.5 gray)
        for rep = 1:maxRep
            hdr_image = img(:,:,((3*rep)-2):(3*rep));
            hdr_image = setExposure(hdr_image,[-log(.5)/median(hdr_image(:)) 1]);
            edges(rep,1) = getElder_edges_KE(uint8(255*hdr_image),nz);
        end
        
        % Reduce file size by saving only the "edge" pixels
        edges = compressElderResults(edges);
        
        % Save results
        save(['SYNS/' int2str(folderNo) '/projected_views/edges_ElderZucker' int2str(nz) '_' int2str(i) '.mat'],'edges');
        clear edges
        
    end
    
end
