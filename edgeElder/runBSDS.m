function runBSDS()

impath = 'C:/Users/kehinger/Desktop/BSDS300/images/test/';
savepath = 'C:/Users/kehinger/Desktop/BSDS300/alg/';

edgewidth = 20;		
edgelplotflag = 0;      % View edgels on image
imageplotflag = 0;      % Don't view imput image
sepplotflag = 0;        % View edgels separately
maxscale = 6;
        
savefilenames = {'edge.mat','blur.mat','dark.mat','light.mat',...
 'g1mag.mat','g1dir.mat','g1scale.mat','g2mag.mat','g2scale.mat',...
 'xzero.mat','yzero.mat','xdark.mat','ydark.mat','xlight.mat','ylight.mat',...
 'Rdark.mat','Gdark.mat','Bdark.mat','Rlight.mat','Glight.mat','Blight.mat',...
 'nxend.mat', 'nyend.mat', 'pxend.mat', 'pyend.mat'};

save_v = zeros(1,length(savefilenames));    % Files to save
view_v = zeros(1,length(savefilenames));   % Files to view
filterdir = 'filters';   
edgetype = '.mat';

imlist = dir([impath '*.jpg']);
for i = 1:length(imlist)
    i
    tic
    imname = imlist(i).name;
    imname = imname(1:end-4);
    img = imread([impath imname '.jpg']);
    
    edgemap = zeros(size(img,1),size(img,2));
    for nz = 1:100
        
        tmp = main_edge([impath imname '.jpg'],maxscale,nz,edgewidth,1,1, ...
            save_v,view_v,filterdir, edgetype, edgelplotflag,imageplotflag,...
            savefilenames,1);
        %{
        disp('Edges found');
        sum(tmp.edge(:))
        %}
        % Should get a subset of previous edges; check for new
        %{
        checkMe = tmp.edge & (edgemap == 0);
        if sum(checkMe(:))>0
            disp([int2str(sum(checkMe(:))) ' new edges at nz = ' int2str(nz) ' in img ' int2str(i)]);
        end
        %}
        edgemap = edgemap + tmp.edge;
        clear tmp
    end
    
    imwrite(edgemap./max(edgemap(:)),[savepath imname '.bmp']);
    toc
end