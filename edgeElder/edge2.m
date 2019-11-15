function varargout = edge(varargin)
% EDGE Application M-file for edge.fig
%    FIG = EDGE launch edge GUI.
%    EDGE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 31-Jul-2002 11:18:28

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    %#######################################
    % Set default condition for convolution:
    %#######################################
    set(handles.radiobutton2, 'Value', 1);
    set(handles.radiobutton3, 'Value', 1);
    
    %##########################################
    % Set default condition for edge file type:
    %##########################################
    set(handles.checkbox42, 'Value', 1);
    
    
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



%##########################################################################
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

	imagefile   = get(handles.edit1,'String');
	scale       = str2num(get(handles.edit2,'String'));
	noise       = str2num(get(handles.edit3,'String'));
	edgewidth   = str2num(get(handles.edit4,'String'));
	con_scale   = get(handles.radiobutton1,'Value');
	con_grad    = get(handles.radiobutton3,'Value');
	outputdir   = get(handles.edit8,'String');
	filterdir   = get(handles.edit12,'String');
    edgelplotflag = get(handles.checkbox46,'Value');
    imageplotflag = get(handles.checkbox49,'Value');
    subpixelflag  = get(handles.checkbox48,'Value');
	
	for i=1:11
        startn  = 20;
        handle_save = strcat('checkbox',num2str(startn+2*(i-1)));
        handle_view = strcat('checkbox',num2str(startn+2*i-1));
        hs = findobj('Tag',handle_save);
        hv = findobj('Tag',handle_view);
        save_v(i) = get(hs,'Value');
        view_v(i) = get(hv,'Value');
	end;
	
    %##################################
    % Check that input files specified:
    %##################################
    if isempty(imagefile)
        mess1   = 'Input image file must be specified. Click ok to re-try.';
        uiwait(msgbox(mess1,'Missing Image Input','warn','modal'));
        return;
    end;
    if isempty(filterdir)
        mess1   = 'Filter directory must be specified. Click ok to re-try.';
        uiwait(msgbox(mess1,'Missing Filter Directory','warn','modal'));
        return;
    end;
    if (isempty(scale) | isempty(noise) | isempty(edgewidth))
        mess1   = ['Scale, noise and edgewidth values must be specified. ',...
                    'Click ok to re-try.'];
        uiwait(msgbox(mess1,'Missing Parameter Values','warn','modal'));
        return;
    end;
    
    %################################################################
	% Check that input filter directory contains the necessary files:
    %################################################################
	ffiles = dir(filterdir);
	fnames = {ffiles(:).name};
	required_filters = {'g1x05.ascii','g1x1.ascii','g1y05.ascii','g1y1.ascii',...
                        'g2x05.ascii','g2x1.ascii','g2y05.ascii','g2y1.ascii',...
                        'gx05.ascii','gx1.ascii','gy05.ascii','gy1.ascii'};
	
    
    memflag = ismember(required_filters,fnames);
    memind  = find(~memflag);
	mess1   = 'Warning: Filter directory is missing the following filter(s):';
    mess3   = [];
    for ind = 1:1:length(memind)
        mess3   = [mess3,', ',required_filters{memind(ind)}];
    end;
	filtmessage = [mess1,mess3]; 
	if ~isempty(memind)    
        uiwait(msgbox(filtmessage,'Filter Directory Warning','warn','modal'));
	end;
    
    %###########################
	% Obtain edge map file type:
    %###########################
    if get(handles.checkbox42,'Value')
        edgetype = '.mat';
    elseif get(handles.checkbox43,'Value')
        edgetype = '.jpg'; 
    elseif get(handles.checkbox44,'Value')
        edgetype = '.tif'; 
    else
        edgetype = '.png'; 
    end;

    %##############################################################
    % Check that output directory name is valid if saving any maps:
    %##############################################################
    [savefilenames] = check_save_requirements(save_v,outputdir,edgetype);

    %#########################
    % Call main_edge function:
    %#########################
	OD = main_edge(imagefile,scale,noise,edgewidth,...
                con_scale,con_grad,save_v,view_v,filterdir,...
                edgetype,edgelplotflag,imageplotflag,savefilenames,...
                subpixelflag);

    assignin('base','Edge_Data',OD);

return;

%##########################################################################
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% 	helpw;  % This is used to call helpw.fig
        uiwait(msgbox('No help files exist.','Help Notice','modal'));
return;

%##########################################################################
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
	delete(handles.figure1);
return;

%##########################################################################
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
	[fname,pname] = uigetfile('*.*');
    if ~(isequal(fname,0) | isequal(pname,0))
    	filename = strcat(pname,fname);
	    set(handles.edit1, 'String', filename);
    end;
return;

%##########################################################################
function varargout = edit1_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = pushbutton14_Callback(h, eventdata, handles, varargin)
  
    sdir = uigetdir('Select Filter Directory',pwd);

    set(handles.edit12,'String',sdir);
    
return;

%##########################################################################
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
	set(handles.edit1,'String','');
	set(handles.edit2,'String','');
	set(handles.edit3,'String','');
	set(handles.edit4,'String','');
	set(handles.edit8,'String','');
	set(handles.edit12,'String','');
	set(handles.checkbox46,'Value',0);
	set(handles.checkbox48,'Value',0);
	set(handles.checkbox49,'Value',0);
return;

%##########################################################################
function varargout = edit2_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = edit3_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = edit4_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = pushbutton6_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.radiobutton1, 'Value');
   set(handles.radiobutton2,'Value',~v1);
return;

%##########################################################################
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.radiobutton2, 'Value');
   set(handles.radiobutton1,'Value',~v1);
return;

%##########################################################################
function varargout = pushbutton1_ButtonDownFcn(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = pushbutton8_Callback(h, eventdata, handles, varargin)

    sdir = uigetdir('Select Output Directory',pwd);

    set(handles.edit8,'String',sdir);


return;

%##########################################################################
function varargout = edit8_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = edit9_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = edit10_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = edit11_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox20_Callback(h, eventdata, handles, varargin)
return;
 
%##########################################################################
function varargout = checkbox21_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox22_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox23_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox24_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox25_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox26_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox27_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox28_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox29_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox30_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox31_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox32_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox33_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox34_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox35_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox36_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox37_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox38_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox38, 'Value');
    if v1
        set(handles.checkbox48,'Value',1);
    end;
return;

%##########################################################################
function varargout = checkbox39_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox39, 'Value');
    if v1
        set(handles.checkbox48,'Value',1);
    end;
return;

%##########################################################################
function varargout = checkbox40_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox40, 'Value');
    if v1
        set(handles.checkbox48,'Value',1);
    end;
return;

%##########################################################################
function varargout = checkbox41_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox41, 'Value');
    if v1
        set(handles.checkbox48,'Value',1);
    end;
return;

%##########################################################################
function varargout = radiobutton3_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.radiobutton3, 'Value');
   set(handles.radiobutton4,'Value',~v1);
return;

%##########################################################################
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.radiobutton4, 'Value');
   set(handles.radiobutton3,'Value',~v1);
return;

%##########################################################################
function varargout = edit12_Callback(h, eventdata, handles, varargin)
return;

%##########################################################################
function varargout = checkbox42_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.checkbox42, 'Value');
    if v1
        set(handles.checkbox43,'Value',0);
        set(handles.checkbox44,'Value',0);
        set(handles.checkbox45,'Value',0);
    else
        set(handles.checkbox43,'Value',1);
        set(handles.checkbox44,'Value',0);
        set(handles.checkbox45,'Value',0);
    end;
return;

%##########################################################################
function varargout = checkbox43_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.checkbox43, 'Value');
    if v1
        set(handles.checkbox42,'Value',0);
        set(handles.checkbox44,'Value',0);
        set(handles.checkbox45,'Value',0);
    else
        set(handles.checkbox42,'Value',1);
        set(handles.checkbox44,'Value',0);
        set(handles.checkbox45,'Value',0);
    end;
return;

%##########################################################################
function varargout = checkbox44_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.checkbox44, 'Value');
    if v1
        set(handles.checkbox42,'Value',0);
        set(handles.checkbox43,'Value',0);
        set(handles.checkbox45,'Value',0);
    else
        set(handles.checkbox42,'Value',1);
        set(handles.checkbox43,'Value',0);
        set(handles.checkbox45,'Value',0);
    end;
return;

%##########################################################################
function varargout = checkbox45_Callback(h, eventdata, handles, varargin)
   v1 = get(handles.checkbox45, 'Value');
    if v1
        set(handles.checkbox42,'Value',0);
        set(handles.checkbox43,'Value',0);
        set(handles.checkbox44,'Value',0);
    else
        set(handles.checkbox42,'Value',1);
        set(handles.checkbox43,'Value',0);
        set(handles.checkbox44,'Value',0);
    end;
return;

%##########################################################################
function varargout = checkbox46_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox46, 'Value');
    if v1
        set(handles.checkbox48,'Value',1);
    end;
return;

%##########################################################################
function varargout = checkbox48_Callback(h, eventdata, handles, varargin)
    v1 = get(handles.checkbox48, 'Value');
    if ~v1
        set(handles.checkbox46,'Value',0);
        set(handles.checkbox38,'Value',0);
        set(handles.checkbox39,'Value',0);
        set(handles.checkbox40,'Value',0);
        set(handles.checkbox41,'Value',0);
    end;
return;

%##########################################################################
function varargout = checkbox49_Callback(h, eventdata, handles, varargin)
return;


%###########################################################################
function[sfilenames] = check_save_requirements(save_flags,outputdir,etype);

    name1       = ['edge',etype];
    sfilenames  = {name1,'blur.mat','dark.mat','light.mat',...
                   'g1mag.mat','g1dir.mat','g1scale.mat','g2mag.mat','g2scale.mat',...
                   'xzero1.mat','yzero1.mat','xzero2.mat','yzero2.mat'};
    
	if any(save_flags)
        
        % Output directory does not exist:
        if ((exist(outputdir) ~= 7) & (length(outputdir) ~= 0))
            mkdir(outputdir);        % Create this directory    
            mess1   = 'Requested output directory has been created.';
            uiwait(msgbox(mess1,'Output Directory Creation','none','modal'));
        end;
        
        if (length(outputdir) == 0)
            outputdir   = 'tmp_output';
            opmessage   = ['Output directory name not specified.  ',...
                           'Files will be saved in directory "tmp_output".'];
            uiwait(msgbox(opmessage,'Output Directory Warning','warn','modal'));
            if (exist(outputdir) ~= 7)
                mkdir(outputdir);
            end;
        end;
            
        filestosave = sfilenames(find(save_flags));
        
        x   = dir(outputdir);
        xx  = {x.name}';
        commonfiles = intersect(filestosave,xx);
    
        if ~isempty(commonfiles)
            maxnumx = 0;
            for k = 1:1:length(xx)
                y = xx{k};
                numx = 0;
                if (length(y)>4)
                    if (y(end-4)=='x')
                        numx = sum(y=='x');
                    end;
                end;
                if (numx > maxnumx)
                    maxnumx = numx;
                end;
            end;
            assignin('base','maxnumx',maxnumx+1);
            assignin('base','commonfiles',commonfiles);
  
            uiwait(outputdir_erase_info);  % Call GUI
        
            changename_flag = evalin('base','changename_flag;');
        
            if changename_flag
               repx = repmat('x',1,maxnumx+1);
               for k = 1:1:length(sfilenames)
                   s    = sfilenames{k};
                   s    = [s(1:end-4),repx,s(end-3:end)];
                   sfilenames{k} = s;
               end;
            end;
           
            evalin('base','clear changename_flag maxnumx commonfiles;');    
        end;    
    end;
    
    for k = 1:1:length(sfilenames)
        s   = sfilenames{k};
        s   = [outputdir,'/',s];
        sfilenames{k} = s;
    end;
    
return;

% %###########################################################################
