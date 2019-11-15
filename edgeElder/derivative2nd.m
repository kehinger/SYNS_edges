%############################################################
% 
% derivative2nd(g1dir,maxscale,noise,gauss_a,conv_type,fpath)
%
%############################################################

function[g2mag,g2sc,g2all] = derivative2nd(g1dir,maxscale,noise,...
                                 gauss_a,conv_type,fpath);

fm2 = '.ascii';

%############
% Initialise:
%############
g2mag   = 0;
g2sc    = 0;
nrows   = size(g1dir,1);
g2all   = zeros(maxscale*nrows,size(g1dir,2));

for scale = 1:1:maxscale
    
    %#############################
    % Set scale value for filters:
    %#############################
    if (scale == 1)
        g2scaleval = '05';
    else
        g2scaleval = '1';
    end;
    
    mimg    = gauss_a(:,:,scale);
    
    kern1   = load(strcat(fpath,'\gx',g2scaleval,fm2));
    rc1     = convolve_2(mimg,kern1,conv_type);
    kern2   = load(strcat(fpath,'\g2y',g2scaleval,fm2));
    rc2     = convolve_2(rc1,kern2,conv_type);
    
    kern3   = load(strcat(fpath,'\gy',g2scaleval,fm2));
    rc3     = convolve_2(mimg,kern3,conv_type);
    kern4   = load(strcat(fpath,'\g2x',g2scaleval,fm2));
    rc4     = convolve_2(rc3,kern4,conv_type);
    
    kern5   = load(strcat(fpath,'\g1x',g2scaleval,fm2));
    rc5     = convolve_2(mimg,kern5,conv_type);
    kern6   = load(strcat(fpath,'\g1y',g2scaleval,fm2));
    rc6     = convolve_2(rc5,kern6,conv_type);
         
    %#######################################
    % Calculate the 2nd Gaussian derivative:
    %#######################################
    g2      = g2steer(rc4,rc2,rc6,g1dir);
     
    g2all((scale-1)*nrows+1:scale*nrows,:) = g2; % For subpixel loc'n.
        
    %##############################################################
    % Augment multi-scale Gaussian directional 2nd derivative maps:
    %##############################################################
    [g2mag,g2sc] = g2scale(g2mag,g2,g2sc,scale,noise,0);
    
end;

return;