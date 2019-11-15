% compressElderResults
%
% Save Elder edge detection results in a more compressed format. The
% original code saves all the results as image maps which are zero
% everywhere an edge was not detected. This format saves only the "edge"
% pixels and discards all the zeros.

function outStruct = compressElderResults(edgeStruct)

for i = 1:size(edgeStruct,1)
    old = edgeStruct(i,1);
    new = [];
    h = size(old.edge,1); % height of image
    sc = size(old.g2_all,1)/h; % number of scales
    ind = find(old.edge > 0);
    new.edge = ind;
    new.blur = old.blur(ind);
    new.dark = old.dark(ind);
    new.light = old.light(ind);
    new.g1mag = old.g1mag(ind);
    new.g1dir = old.g1dir(ind);
    new.g1scale = old.g1scale(ind);
    new.g2mag = old.g2mag(ind);
    new.g2scale = old.g2scale(ind);
    % g2_all is results from each scale stacked in dim 1
    new.g2_all = zeros(length(ind),sc);
    for s = 1:sc
        panel = old.g2_all(((s-1)*h)+1:(s*h),:);
        new.g2_all(:,s) = panel(ind);
    end
    new.xzero = old.xzero(ind);
    new.yzero = old.yzero(ind);
    if isfield(old,'Rdark')
        new.Rdark = old.Rdark(ind);
    end
    if isfield(old,'Gdark')
        new.Gdark = old.Gdark(ind);
    end
    if isfield(old,'Bdark')
        new.Bdark = old.Bdark(ind);
    end
    if isfield(old,'Rlight')
        new.Rlight = old.Rlight(ind);
    end
    if isfield(old,'Glight')
        new.Glight = old.Glight(ind);
    end
    if isfield(old,'Blight')
        new.Blight = old.Blight(ind);
    end
    new.nxend = old.nxend(ind);
    new.nyend = old.nyend(ind);
    new.pxend = old.pxend(ind);
    new.pyend = old.pyend(ind);
    outStruct(i,1) = new;
end