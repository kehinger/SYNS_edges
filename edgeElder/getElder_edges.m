% Functional form of James' edge extraction algorithm
%
% function getElder_edges(ImgName,noise)
%
% ImgName: name of the file to process
% nz: noise level in the image -also controls detail level-
%
% All output is saved to the current working directory
%
%
% K.E. notes:
% - Moved output files to a subfolder (automatically created)
% - This code assumes 3 color channels
% - Has image boundary artifacts (leave some padding around input patch)
% - Increasing nz or increasing edgewidth reduces detection of smaller /
%   blurrier edges (what is the difference between these?)
% - Maxscale works like you'd expect, it adds edges at increasing spatial
%   scale (blurrier edges)
% - Minscale isn't used


function OD = getElder_edges(ImgName,nz,minscale,maxscale)

% Extract the edges from the image
edgewidth = 20;		
edgelplotflag = 1;      % View edgels on image
imageplotflag = 1;      % Don't view imput image
sepplotflag = 1;        % View edgels separately               
        
savefilenames = {'edge.mat','blur.mat','dark.mat','light.mat',...
 'g1mag.mat','g1dir.mat','g1scale.mat','g2mag.mat','g2scale.mat',...
 'xzero.mat','yzero.mat','xdark.mat','ydark.mat','xlight.mat','ylight.mat',...
 'Rdark.mat','Gdark.mat','Bdark.mat','Rlight.mat','Glight.mat','Blight.mat',...
 'nxend.mat', 'nyend.mat', 'pxend.mat', 'pyend.mat'};

% Save results in a subfolder named after the input image - KE, Sept 2016
k = strfind(ImgName,'/');
if isempty(k)
    n = ImgName;
else
    n = ImgName(k(end)+1:end);
end
n = n(1:end-4);
savefilenames = strcat([n '/'],savefilenames);
mkdir(n);

save_v = ones(1,length(savefilenames));    % Files to save
view_v = zeros(1,length(savefilenames));   % Files to view
view_v(1) = 1;
filterdir = 'filters/';   
edgetype = '.mat';

% Run the main edge extraction program
OD = main_edge(ImgName, maxscale,nz,edgewidth,1,1, ...
    save_v,view_v,filterdir, edgetype, edgelplotflag,imageplotflag,...
    savefilenames,1);   

%close all;

return;