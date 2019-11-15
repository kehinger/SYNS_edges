function OD = getElder_edges_KE(img,varargin)

if isempty(varargin)
    nz = 10.24;
else
    nz = varargin{1};
end

edgewidth = 20;
minscale = 1;
maxscale = 6;

savefilenames = {'edge.mat','blur.mat','dark.mat','light.mat',...
 'g1mag.mat','g1dir.mat','g1scale.mat','g2mag.mat','g2scale.mat',...
 'xzero.mat','yzero.mat','xdark.mat','ydark.mat','xlight.mat','ylight.mat',...
 'Rdark.mat','Gdark.mat','Bdark.mat','Rlight.mat','Glight.mat','Blight.mat',...
 'nxend.mat', 'nyend.mat', 'pxend.mat', 'pyend.mat'};

save_v = zeros(1,length(savefilenames));    % Files to save
view_v = zeros(1,length(savefilenames));   % Files to view
filterdir = 'filters/';   
edgetype = '.mat';

% Run the main edge extraction program
OD = main_edge(img, maxscale, nz, edgewidth, 1, 1, ...
    save_v, view_v, filterdir, edgetype, 0, 0, savefilenames,1);

end