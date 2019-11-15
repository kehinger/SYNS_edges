%############################################################
%
% main_edge.m:  Main function called from edge GUI interface.
%
%############################################################

function[OutputData] = main_edge(imagefile,scale,noise,edgew,...
                        conss,congrd,save_v,view_v,filtpath,...
                        edgetype,edgelplotflag,imgplotflag,savenames,...
                        subpixelflag)

%#############
% Input Image:
%#############

% Modified this to accept either an image matrix or a filename as the first
% input to the function. - KE, Sept 2016
% Added a "noround" version of the C code that does not round the
% light/dark maps to integers (but it does still truncate values above
% 255). - KE, Jan 2017

if ~isa(imagefile,'numeric')
    ind1    = strfind(imagefile,'.');
    ftype   = imagefile(ind1+1:end);
    
    if (strcmp(ftype,'ascii') == 1)
        img = dlmread(imagefile,' ');
    else
        img = imread(imagefile);
        
    end
    imagetitle = imagefile;
else
    img = imagefile;
    imagetitle = 'IMAGE';
end

if (size(img,3)~=1)
    imgR = img(:, :, 1);
    imgG = img(:, :, 2);
    imgB = img(:, :, 3);
    img = rgb2gray(img);
else
    imgR = img;
    imgG = img;
    imgB = img;
end
img=double(img);


path        = pwd;
dirind      = strfind(path,'\');
datapath    = strcat(path(1:dirind(end)),'\data\');

%########################
% Compose blurred images:
%########################
gauss_imgs  = scalespace(img,scale,conss);

%########################
% Calculate gradient map:
%########################
[g1_mag,g1_dir,g1_sc] = gradient(scale,noise,gauss_imgs,congrd,filtpath);

%##############################
% Calculate 2nd derivative map:
%##############################
[g2_mag,g2_sc,g2_all] = derivative2nd(g1_dir,scale,noise,gauss_imgs,congrd,filtpath);

fprintf('\n**** Starting Localization Routines ****\n');
tic;
%####################
% Calculate edge map:
%####################

clear find_edges_fc
[edge_map,blur_map,dark_map,light_map,xzero1_map,yzero1_map,xzero2_map,yzero2_map, ...
  Rdark_map,Gdark_map,Bdark_map,Rlight_map,Glight_map,Blight_map, nxend_map, nyend_map, pxend_map, pyend_map] = ...
    find_edges_fc_noround(img,g1_mag,g1_dir,g1_sc,g2_mag,g2_sc,g2_all,noise,edgew,subpixelflag,scale,double(imgR),double(imgG),double(imgB));
clear find_edges_fc     

% Branka's code leaves dark and light as double which need rounding:
% Nope, light and dark are integers if you use find_edges_fc. - KE
% dark_map    = round(dark_map);
% light_map    = round(light_map);
xzero_map=xzero1_map;
yzero_map=yzero1_map;

toc3 = toc;
fprintf('\nTime taken by CVPR Algorithm: %f\n',toc3);

%dark_map   = round(dark_map);
%light_map  = round(light_map);

% NOTE: These maps are already mostly integers, but may include some small
% double values.
Rdark_map=round(Rdark_map);
Gdark_map=round(Gdark_map);
Bdark_map=round(Bdark_map);
Rlight_map=round(Rlight_map);
Glight_map=round(Glight_map);
Blight_map=round(Blight_map);


%#######################
% Output data in struct:
%#######################
OutputData = struct('edge',edge_map,'blur',blur_map,...
    'dark',dark_map,'light',light_map,...
    'g1mag',g1_mag,'g1dir',g1_dir,'g1scale',g1_sc,...
    'g2mag',g2_mag,'g2scale',g2_sc,'g2_all',g2_all,...
    'xzero',xzero_map,'yzero',yzero_map,...
    'Rdark',Rdark_map,'Gdark',Gdark_map,'Bdark',Bdark_map,...
    'Rlight',Rlight_map,'Glight',Glight_map,'Blight',Blight_map, 'nxend',...
    nxend_map, 'nyend', nyend_map, 'pxend', pxend_map, 'pyend', pyend_map);

fprintf('\n**** Finished Localization Routines ****\n');

%#####################
% View original image:
%#####################
if imgplotflag
    figure;
    colormap(gray);
    axis image;
    imagesc(img);
    axis image;
end;

%###############
% View edge map:
%###############
if view_v(1)
    figure;
    colormap(gray);
    axis image;
    imagesc(edge_map);
    axis image;
end;

%###############
% Save edge map:
%###############
if save_v(1)
    if strcmp(edgetype,'.mat')
        save(savenames{1},'edge_map');
    else
		imwrite(uint8(edge_map),savenames{1});
    end;
end;
        
%#########################
% View edgel map on image:
%#########################
if edgelplotflag
    newarrcolor = [1 0 1];
    legtext1    = 'Centred edgels';
    legtext2    = 'Subpixel method 1';
    legtext3    = 'Subpixel method 2';
    plot_edgelmap(img,edge_map,xzero1_map,yzero1_map,g1_dir,...
                [1 0 1],legtext1,legtext2,imagetitle);
    plot_edgelmap(img,edge_map,xzero2_map,yzero2_map,g1_dir,...
                [1 0 0],legtext1,legtext3,imagetitle);
end;

%########################
% View maps as requested:
%########################
if view_v(2)
    figure; imagesc(blur_map); axis image; title('Blur');
end;
if view_v(3)
    figure; imagesc(dark_map); axis image; title('Dark');
end;
if view_v(4)
    figure; imagesc(light_map); axis image; title('Light');
end;
if view_v(5)
    figure; imagesc(g1_mag); axis image; title('Gradient Magnitude');
end;
if view_v(6)
    figure; imagesc(g1_dir); axis image; title('Gradient Direction');
end;
if view_v(7)
    figure; imagesc(g1_sc); axis image; title('Gradient Scale');
end;
if view_v(8)
    figure; imagesc(g2_mag); axis image; title('2nd Derivative Magnitude');
end;
if view_v(9)
    figure; imagesc(g2_sc); axis image; title('2nd Derivative Scale');
end;
if view_v(10)
    figure; imagesc(xzero1_map); axis image; title('X zero-crossing offset - method 1');
    figure; imagesc(xzero2_map); axis image; title('X zero-crossing offset - method 2');
end;
if view_v(11)
    figure; imagesc(yzero_map); axis image; title('Y zero-crossing offset - method 1');
    figure; imagesc(yzero_map); axis image; title('Y zero-crossing offset - method 2');
end;

%########################
% Save maps as requested:
%########################
if save_v(1)
    save(savenames{1}, 'edge_map');
end
if save_v(2)
    save(savenames{2},'blur_map');
end;
if save_v(3)
    save(savenames{3},'dark_map');
end;
if save_v(4)
    save(savenames{4},'light_map');
end;
if save_v(5)
    save(savenames{5},'g1_mag');
end;
if save_v(6)
    save(savenames{6},'g1_dir');
end;
if save_v(7)
    save(savenames{7},'g1_sc');
end;
if save_v(8)
    save(savenames{8},'g2_mag');
end;
if save_v(9)
    save(savenames{9},'g2_sc');
end;
if save_v(10)
    save(savenames{10},'xzero_map');
end;
if save_v(11)
    save(savenames{11},'yzero_map');
end;
if save_v(16)
    save(savenames{16},'Rdark_map');
    save(savenames{17},'Gdark_map');
    save(savenames{18},'Bdark_map');
    save(savenames{19},'Rlight_map');
    save(savenames{20},'Glight_map');
    save(savenames{21},'Blight_map');
end
if save_v(22)
    save(savenames{22}, 'nxend_map');   
    save(savenames{23}, 'nyend_map'); 
    save(savenames{24}, 'pxend_map'); 
    save(savenames{25}, 'pyend_map'); 
end

return;
