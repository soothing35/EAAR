
function varargout = spiral_picker_GUI(varargin)
% SPIRAL_PICKER_GUI M-file for spiral_picker_GUI.fig
%      SPIRAL_PICKER_GUI, by itself, creates a new SPIRAL_PICKER_GUI or raises the existing
%      singleton*.
%
%      H = SPIRAL_PICKER_GUI returns the handle to a new SPIRAL_PICKER_GUI or the handle to
%      the existing singleton*.
%
%      SPIRAL_PICKER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIRAL_PICKER_GUI.M with the given input arguments.
%
%      SPIRAL_PICKER_GUI('Property','Value',...) creates a new SPIRAL_PICKER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spiral_picker_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spiral_picker_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spiral_picker_GUI

% Last Modified by GUIDE v2.5 11-Jan-2013 11:30:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spiral_picker_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @spiral_picker_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spiral_picker_GUI is made visible.
function spiral_picker_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spiral_picker_GUI (see VARARGIN)

 
%%%%% always show this main window in the center
%pixels
set( gcf, 'Units', 'pixels' );


%get your display size
screenSize = get(0, 'ScreenSize');

%calculate the center of the display
position = get( gcf, 'Position' );
position(1) = (screenSize(3)-position(3))/2;
position(2) = (screenSize(4)-position(4))/2-30;
  
%center the window
set( gcf, 'Position', position );

% set the button figures
set(handles.pushbutton_nextName,'cdata',arrowbutton);
set(handles.pushbutton_formerName,'cdata',inarrowbutton);

% initiate some data of handles and some signs
[Newhandles]= initiate_handles(handles);
handles= Newhandles;


%%% read a parameter file to get the parameters and image directory
% Find our local directory
pa=fileparts(which('spiral_picker_GUI'));
% Retrieve parameters from BoxPara.mat in the local directory
if exist([pa '/SpiralPara.mat'],'file')
    load([pa '/SpiralPara.mat']);
    if isempty(scale)||isempty(currentfile_path)
        scale=1;currentfile_path=pa;
    end
else
    scale = 1;currentfile_path =pa;
    
end
if ~isempty(currentfile_path)
    handles.currentfile_path = currentfile_path;
else
    handles.currentfile_path = '';
end


handles.procPath ='';
handles.InfoPath ='';
handles.resultPath = '';
handles.filename = '';
handles.fullPath ='';


handles.pickednum =0;
handles.scale = 1;

%set(gcf,'windowbuttonmotionfcn','if ~isempty(overobj(''handles.axes1'')) disp(''hello''); end');

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spiral_picker_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spiral_picker_GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_formerName.
function pushbutton_formerName_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_formerName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)






% --- Executes on button press in pushbutton_loadfile.
function pushbutton_loadfile_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.currentfile_path)&& exist(handles.currentfile_path,'dir')
    [FileName,PathName] = uigetfile( fullfile(handles.currentfile_path,'*.mrc;*.tif'),'Load image File');
else
    [FileName,PathName] = uigetfile({'*.mrc','*.tif'},'Load image File');
end
if (FileName==0)
    return;
elseif (~strcmp(FileName(end-2:end),'mrc'))&& (~strcmp(FileName(end-2:end),'tif'))  % judge whether choose image file
    msgbox('Wrong image file type!','Error','error');
    return;
end

a=strfind(PathName,'/'); str= PathName(1: a(size(a,2)-1));
%a=strfind(PathName,'\'); str= PathName(1: a(size(a,2)-1));

handles.currentfile_path = str;
handles.procPath =[handles.currentfile_path 'Image/'];
%handles.procPath =[handles.currentfile_path 'Image\'];

%     handles.InfoPath =[handles.currentfile_path 'Info/'];
handles.resultPath = [handles.currentfile_path 'Picking-result/'];
%handles.resultPath = [handles.currentfile_path 'Picking-result\'];

if ~exist([handles.resultPath],'dir')
    mkdir([handles.resultPath] )
end

    
[~, b, Ext] = fileparts(FileName);
handles.baseName = b;
    
handles.fullPath = [handles.procPath FileName];  
     
 
 
handles.filename = FileName;
imagefile_ext = Ext(strfind(Ext,'.')+1:end);


% generate mat file name for saving
% find the image name
a=strfind(handles.filename,'.'); 
str= handles.filename(1: a-1);

% find the points files' index number
dirData = dir([handles.resultPath, '*.mat']);    %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
% if ~isempty(fileList)
%     fileList = cellfun(@(x) fullfile(handles.resultPath,x),...  %# Prepend path to files
%                    fileList,'UniformOutput',false);
% end

handles.nextfile_Num = length(fileList)+1;
newfile_Num = num2str(handles.nextfile_Num); 
handles.save_matfile = [handles.resultPath str '_' newfile_Num '.mat'] ;


% %%%% if reload a new file, we need to clear some parameters
% initiate some data of handles and some signs
[Newhandles]= initiate_handles(handles);
handles= Newhandles;

[Newhandles]=show_image(hObject,imagefile_ext,handles);
handles=Newhandles;

 

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in pushbutton_nextName.
function pushbutton_nextName_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_nextName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




 



% --- Executes on button press in pushbutton_contrast.
function pushbutton_contrast_Callback(~, ~, handles)
% hObject    handle to pushbutton_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.fileLoaded   
    imcontrast(handles.axes1);
end


% --- Executes on selection change in popupmenu_filters_lowpass.
function popupmenu_filters_lowpass_Callback(hObject, ~, handles)
% hObject    handle to popupmenu_filters_lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_filters_lowpass contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_filters_lowpass
if  handles.fileLoaded          % if already loaded file

    contents = cellstr(get(hObject,'String')); 
    filter_name =  contents{get(hObject,'Value')};

    button_state = get(handles.togglebutton_invertImage,'Value');
    if button_state == get(handles.togglebutton_invertImage,'Max')
        % Toggle button is pressed, take appropriate action
        a= size(handles.originalim);
        white_mat = uint8(255.*ones(a));
        handles.currentim = white_mat-handles.originalim;
    else  
        handles.currentim = handles.originalim;
 
    end
 

    RUN = 1;
    if (strcmp(filter_name,'none'))
        if (handles.scale ==1)
            set(handles.ih,'Cdata',handles.currentim);
            handles.Grayim = handles.currentim;
        else
            %Scaleim = uint8(BinImage(handles.originalim,handles.scale,'mean'));
            Scaleim = uint8(imresize(handles.currentim,1/handles.scale ));
            handles.Grayim = Scaleim;     
        end
    
     elseif (strcmp(filter_name,'Gaussian'))     % Gaussian filter
        while (RUN==1)
%             prompt = {'Enter standard deviation sigma (positive):'};
            prompt = {'Enter corner frequency (positive):'};

            dlg_title = 'Enter Gaussian Parameters:';
            num_lines = 1;
%            def = {'0.2'};
            def = {'0.1'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if (isempty(answer))
                set(handles.popupmenu_filters_lowpass,'Value',1);
                return;
            end
            handles.hsize = 10;
%            handles.sigma = str2double(answer{1});

            handles.sigma = sqrt(log(2))/(2*pi*str2double(answer{1})); 
            if ((str2double(answer{1})>=0) && (str2double(answer{1})<=15))
                RUN = 0;
            end
        end

        H = fspecial('gaussian', handles.hsize, handles.sigma);
        w = waitbar(0, 'Gaussian filtering ... Please wait ...');
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
        gaussian_im=handles.Grayim;
        close(w);
    elseif (strcmp(filter_name,'Average'))       % average filter
        while (RUN==1)
            prompt = {'Enter Average filter size (positive):'};
            dlg_title = 'Enter Average Parameters:';
            num_lines = 1;
            def = {'3'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if (isempty(answer))
                set(handles.popupmenu_filters_lowpass,'Value',1);
                return;
            end
            handles.hsize = str2double(answer{1});
            if (str2double(answer{1})>0)
                RUN = 0;
            end
        end        
        H = fspecial('average', handles.hsize);
        w = waitbar(0, 'Average filtering ... Please wait ...');
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
        close(w);       
    end
    
    Newhandles = redraw_image_inMainwin(0,0,handles);  % add filters but no points
    handles  = Newhandles;

    % Update handles structure
    guidata(hObject, handles);

else
    msgbox('No image file has been loaded!','Error','error');
    set(hObject,'Value',1);
    
end

 


% --- Executes during object creation, after setting all properties.
function popupmenu_filters_lowpass_CreateFcn(hObject, ~, handles)
% hObject    handle to popupmenu_filters_lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_filters_highpass.
function popupmenu_filters_highpass_Callback(hObject, ~, handles)
% hObject    handle to popupmenu_filters_highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_filters_highpass contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_filters_highpass
if  handles.fileLoaded

    contents = cellstr(get(hObject,'String')); 
    filter_name =  contents{get(hObject,'Value')};

    button_state = get(handles.togglebutton_invertImage,'Value');
    if button_state == get(handles.togglebutton_invertImage,'Max')
        % Toggle button is pressed, take appropriate action
        a= size(handles.originalim);
        white_mat = uint8(255.*ones(a));
        handles.currentim = white_mat-handles.originalim;
    else  
        handles.currentim = handles.originalim;
 
    end
 

    
    RUN = 1;
    if (strcmp(filter_name,'none'))
        if (handles.scale ==1)
            set(handles.ih,'Cdata',handles.currentim);
            handles.Grayim = handles.currentim;
        else
            %Scaleim = uint8(BinImage(handles.originalim,handles.scale,'mean'));
            Scaleim = uint8(imresize(handles.currentim,1/handles.scale ));
            handles.Grayim = Scaleim;     
        end
    
    elseif (strcmp(filter_name,'Wiener'))
        while (RUN==1)
            prompt = {'Enter Wiener filter size 1 (positive):','Enter Wiener filter size 2 (positive): '};
            dlg_title = 'Enter Wiener Parameters:';
            num_lines = 1;
            def = {'3','3'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if (isempty(answer))
                set(handles.popupmenu_filters_highpass,'Value',1);
                return;
            end
        
            handles.wsize = str2double(answer{1});
            handles.hsize = str2double(answer{2});

            if ((str2double(answer{1})>0) && (str2double(answer{2})>0))
                RUN = 0;
            end
        end
        w = waitbar(0, 'Wiener filtering ... Please wait ...');
        handles.Grayim = wiener2(handles.Grayim,[handles.wsize handles.hsize]);
        close(w);       
    elseif (strcmp(filter_name,'Sharp'))
        while (RUN==1)

            prompt = {'Enter sharp  filter alpha(0-1,positive): '};
            dlg_title = 'Enter sharp Parameters:';
            num_lines = 1;
            def = {'0.2'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if (isempty(answer))
                set(handles.popupmenu_filters_highpass,'Value',1);
                return;
            end
        
            alpha = str2double(answer{1});

            if ((str2double(answer{1})>=0) && (str2double(answer{1})<=1))
                RUN = 0;
            end
        end
        H = fspecial('unsharp');
        w = waitbar(0, 'Unsharp contrast enhancement filter filtering ... Please wait ...');
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
        close(w);       
    end
    
    Newhandles = redraw_image_inMainwin(0,0,handles);   % add filters but no points 
    handles  = Newhandles;
else
    msgbox('No image file has been loaded!','Error','error');
    set(hObject,'Value',1);
end
% Update handles structure
guidata(hObject, handles);


 


% --- Executes on selection change in popupmenu_Scale.
function popupmenu_Scale_Callback(hObject, ~, handles)
% hObject    handle to popupmenu_Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Scale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Scale
str_num = get(hObject,'Value') ;
popStr = get(hObject,'string') ;
p = popStr(str_num);
Newscale= cellfun(@(x)str2double(x), p);

     
if (Newscale~=handles.scale) && handles.fileLoaded
    handles.scale = Newscale;              % actually, the scale is 1/Newscale, we only store the denominator
    
    axes(handles.axes1);
    delete(handles.ih);      % delete the original one
    
    button_state = get(handles.togglebutton_invertImage,'Value');
 
 if button_state == get(handles.togglebutton_invertImage,'Max')
	% Toggle button is pressed, take appropriate action
    a= size(handles.originalim);
    white_mat = uint8(255.*ones(a));
    handles.currentim = white_mat-handles.originalim;
 else  
    handles.currentim = handles.originalim;
 
 end
 
    if (Newscale ==1)
        ih = imshow(handles.currentim,'InitialMagnification',100,'DisplayRange',[min(get(gca,'CLim')) max(get(gca,'CLim'))]);
        handles.Grayim = handles.currentim;

    else
        %Scaleim = uint8(BinImage(handles.originalim,Newscale,'mean'));
        Scaleim = uint8(imresize(handles.currentim,1/Newscale ));

        ih = imshow(Scaleim,'DisplayRange',[min(get(gca,'CLim')) max(get(gca,'CLim'))]);
        handles.Grayim = Scaleim;     
    end    
    
    
    set(handles.editSize, 'String', sprintf('SIZE  : %d x %d', size(handles.Grayim,1),size(handles.Grayim,2)));
     
     panel = findobj(gcf,'Tag','Main_im'); 
     set(handles.axes1,'Parent',panel);
     delete(handles.hSP);
     hSP = imscrollpanel(panel,ih);
     set(hSP,'Units','normalized',...
         'Position',[0 0 1 1]);
     
     hNav = imoverviewpanel(gcf,ih);
     set(hNav,'Units','Normalized','Position',[0.82 0.04 0.19 .19] )

 
     handles.ih=ih;
     handles.hSP=hSP;
     handles.hNav = hNav;
     
     [Newhandles]=redraw_image_inMainwin(1,1,handles);
     handles= Newhandles;
   
     
     % Update handles structure
     guidata(hObject, handles);
     
     % return focus to the main window
     set(hObject, 'Enable', 'off');
     drawnow;
     set(hObject, 'Enable', 'on');
end

 
 





 


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if handles.fileLoaded
    currentfile_path = handles.procPath;   %  save the upper level folder of info and Merge. For example:/Users/yunhui/cryo_em/LiguoData/Data110908
else
    currentfile_path =  handles.currentfile_path;
end

pa=fileparts(which('spiral_picker_GUI'));
scale= handles.scale; 
save([pa '/SpiralPara.mat'],'scale','currentfile_path');

% Construct a questdlg with two options
choice = questdlg('Would you like to quilt?', ...
 'Your choice','OK','No','OK');
% Handle response,if havenot saved, call the "save" function
if strcmp(choice,'OK')
    if (handles.fileLoaded) && (handles.save ==0)
        pushbutton_save_Callback(hObject, [], handles);
    end
    delete(hObject);
end
 

 


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global num;  % points num in one file
global num_all; % all points num in files 


ty=get(gcf,'SelectionType');


switch ty
  case 'normal' % left button was clicked
      % code desired when left button was used
    if (handles.fileLoaded == 1) && (num>=0)
        num = num+1;
        num_all = num_all+1;
        p= round(get(handles.axes1,'CurrentPoint')); %p(1,1) is cols in image matrix(X coordinate)  p(1,2) is rows in image matrix(Y coordinate),matrix origin is on the left top

        api = iptgetapi(handles.hSP);  
        Vpart_pos= api.getVisibleLocation();

        pos = getpixelposition(handles.Main_im);
        scale = handles.scale;
        hightim = size(handles.originalim,1)/scale;
        widthim = size(handles.originalim,2)/scale;
     
        % judge whether the user clicks the image area
        if (p(1,2) < (Vpart_pos(2)+min(pos(4),hightim))) && (p(1,2) >Vpart_pos(2)) &&(p(1,1) < (Vpart_pos(1)+min(pos(3),widthim))) && (p(1,1) >Vpart_pos(1)) 
 
            str_new_x = num2str(p(1,1)*scale);    
            str_new_y = num2str(p(1,2)*scale);
            str_num = num2str(num);
    
            new_string = [str_num ' ' str_new_x ' ' str_new_y];
            
            %show the new clicked cooridinates in the coordinates listbox
            Val= cellstr(get(handles.pixel_Coordinate,'string'));
            NewVal= [Val; {new_string}];
            set(handles.pixel_Coordinate,'string', NewVal);
    
            % save the position of one file
            handles.Points_position(num,1)=p(1,1)*scale;
            handles.Points_position(num,2)=p(1,2)*scale;

            % save all the points positions of all files
            handles.AllPoints_position(num_all,1)=p(1,1)*scale;
            handles.AllPoints_position(num_all,2)=p(1,2)*scale;

            set(handles.Points_num,'string',str_num);
            
            contents = cellstr(get(handles.popupmenu_color,'String'));
            val_color=contents{get(handles.popupmenu_color,'Value')};   % get color
             
            handles.Points_color(num)=val_color(1);
            handles.AllPoints_color(num_all)=val_color(1);

            axes(handles.axes1);

            hold on;
            handles.crossPoints(num)=plot(p(1,1),p(1,2),[val_color(1) '*']);
        end
   
    end
    case 'alt'  % right clicked
        % clear listbox and save the picked points
        set(handles.pixel_Coordinate,'string', '');
          
%         %%%%% change the origin from left top to left bottom. Y = hight of image - Y : 
%         %%%%%% handles.Points_position(:,2) = handles.info(2)-handles.Points_position(:,2);
%     
%         handles.Points_position(:,2) = handles.info(2)-handles.Points_position(:,2);

        Points_position = handles.Points_position;
        Points_color= handles.Points_color;
        a=strfind(handles.filename,'.'); 
        str= handles.filename(1: a-1);


        newfile_Num = num2str(handles.nextfile_Num);
% 
        handles.save_matfile = [handles.resultPath str '_' newfile_Num '.mat'] ;
        
        save(handles.save_matfile, 'Points_position','Points_color');
       
        % clear  handles.Points_position and  handles.Points_color
        handles.Points_position=[];
        handles.Points_color='';
        
        num =0;
        handles.save =0;
        handles.nextfile_Num=handles.nextfile_Num+1;
        newfile_Num = num2str(handles.nextfile_Num);

        handles.save_matfile = [handles.resultPath str '_' newfile_Num '.mat'] ;
        
        set(handles.Points_num, 'String', sprintf('Points Num: %d', num));
        
    case 'extend'    % shift right click : delete the current point
        if num>0
            delete(handles.crossPoints(num));
            handles.crossPoints = handles.crossPoints(1:end-1);
            handles.Points_position = handles.Points_position(1:end-1,:);
            handles.Points_color = handles.Points_color(1:end-1);
            handles.AllPoints_position = handles.AllPoints_position(1:end-1,:);
            handles.AllPoints_color = handles.AllPoints_color(1:end-1);
 
            
        	num = num-1;
            num_all = num_all-1;
        
            set(handles.Points_num, 'String', sprintf('Points Num: %d', num));
            B=get(handles.pixel_Coordinate,'String');
            len=length(B);
        
            if len > 0
                B =  B(1:len-1);
                set(handles.pixel_Coordinate,'String',B);
            end
        end
end
hold off;
% Update handles structure
guidata(hObject, handles);




% --- Hand-written callback 
% --- Used to return 'CData' for the Stop icon on the Record\Stop toggle button
function invarrow = inarrowbutton

invarrow = iconize(imread('inverse_arrow.jpg'));
invarrow(invarrow==255) = .8*255;

 


 % --- Hand-written callback 
% --- Used to return 'CData' for the Play icon on the Play button
function arrow = arrowbutton

arrow = iconize(imread('arrow.jpg'));
arrow(arrow==255) = .8*255;

% --- Hand-written callback 
% --- Used to create icon data from an image, a
function out = iconize(a)

% Find the size of the acquired image and determine how much data will need
% to be lost in order to form a 18x18 icon
[r,c,d] = size(a);
r_skip = ceil(r/9);
c_skip = ceil(c/9);

% Create the 18x18 icon (RGB data)
out =  a(1:r_skip:end,1:c_skip:end,:);


 



 



function [Newhandles]=show_image(~,imagefile_ext,handles)
     


set(handles.editSize, 'Visible', 'on');
    
if strcmp(imagefile_ext,'tif')
    info = imfinfo(handles.fullPath);
    imWidth = info.Width; imHeight=info.Height;
    original =  imread(handles.fullPath);   
    
else   % '.mrc' 'jpg'
    [original s]=ReadEMFile(handles.fullPath);
    imWidth = size(original,1); imHeight = size(original,2);
end
%original = rot90(original);   % note: the origin of the image is in the left bottom. So rotate the matrix 90 degree counter clockwise

set(handles.editSize, 'String', sprintf('SIZE  : %d x %d', round(imWidth/handles.scale), round(imHeight/handles.scale)));
set(handles.textFileName, 'String', handles.filename); 


handles.info = [imWidth imHeight];     % save image's width and height in pixel 
 
Grayim =  uint8(imscale(original));
handles.originalim = Grayim;    % handles.origalim stores the residual image with scale=1
handles.currentim = Grayim;     % handles.currentim stores the current kind of image: vesicle image or residual image 

handles.Grayim = Grayim;        % handles.Grayim sotres the current residual image in main window
handles.fileLoaded = 1;

%set(handles.axes1,'Visible','off'); 
delete (get(handles.axes1,'child'));

panel = findobj(gcf,'Tag','Main_im');
set(handles.axes1,'Parent',panel);
% get the scale parameter
scale =  get(handles.popupmenu_Scale,'Value') ;
% show image according to scale
if (scale ==1)
    ih = imshow(Grayim,'InitialMagnification',100 );
else
    % Scaleim = uint8(BinImage(Grayim,scale,'mean'));
    Scaleim =  imresize(Grayim,1/scale);
    ih = imshow(Scaleim );
    handles.Grayim = Scaleim;
end


handles.ih=ih;
set(handles.ih,'Parent',handles.axes1);

hSP = imscrollpanel(panel,ih);
set(hSP,'Units','normalized','Position',[0 0 1 1]);
    
hNav = imoverviewpanel(gcf,ih);
set(hNav,'Units','Normalized','Position',[0.82 0.04 0.19 .19] )
axes(handles.axes1);

handles.hSP = hSP;
handles.hNav = hNav;
% Newhandles= handles;


Newhandles= handles;  % output parameter         
 

     
function disableButtons(handles)
%change the mouse cursor to an hourglass
set(handles.figure1,'Pointer','watch');

%disable all the buttons so they cannot be pressed

set(handles.pushbutton_loadfile,'Enable','off');
set(handles.pushbutton_formerName,'Enable','off');
set(handles.pushbutton_nextName,'Enable','off');
set(handles.pushbutton_contrast,'Enable','off');

set(handles.popupmenu_filters_lowpass,'Enable','off');
set(handles.popupmenu_filters_highpass,'Enable','off');
set(handles.popupmenu_Scale,'Enable','off');

 

function enableButtons(handles)
%change the mouse cursor to an arrow
set(handles.figure1,'Pointer','arrow');

%enable all the buttons so they can be pressed

set(handles.pushbutton_loadfile,'Enable','on');
set(handles.pushbutton_formerName,'Enable','on');
set(handles.pushbutton_nextName,'Enable','on');
set(handles.pushbutton_contrast,'Enable','on');

set(handles.popupmenu_filters_lowpass,'Enable','on');
set(handles.popupmenu_filters_highpass,'Enable','on');
set(handles.popupmenu_Scale,'Enable','on');





 

 
 
 


 

% --- Executes on key press with focus on pushbutton_contrast and none of its controls.
function pushbutton_contrast_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_contrast (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
k= eventdata.Key; %k is the key that is pressed

if strcmp(k,'return') %if enter was pressed
    pushbutton_contrast_Callback(hObject, [], handles);
end


 
 
 


 

 
 

function  [Newhandles]= initiate_handles(handles)
 
global num;
global num_all;

num_all=0;

num=0;
handles.fileLoaded = 0;   % a sign for loading image file
handles.save = 0;       % a sign for clicked save button
handles.Points_position=[];
handles.Points_color=[];

set(handles.pixel_Coordinate,'string',''); 
set(handles.axes1,'Visible','off');
set(handles.editSize, 'Visible', 'off');
 


set(handles.Points_num, 'String', sprintf('Points Num: %d', num));
set(handles.Points_num, 'Value', num);

Newhandles = handles; 




function h = sfigure(varargin)
    
if nargin>=1
    if ishandle(varargin{1})
        set(0, 'CurrentFigure', varargin{1});
        h=varargin{1};
    else
        h = figure(varargin{:});
    end
else
    h = figure;
end

 


% --- Executes on selection change in popupmenu_color.
function popupmenu_color_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_color


% --- Executes during object creation, after setting all properties.
function popupmenu_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.Points_position,1)~=0 && size(handles.Points_color,1)~=0
    handles.save = 1;
    
    
%     %%%%% change the origin from left top to left bottom. Y = hight of image - Y : 
%     %%%%%% handles.Points_position(:,2) = handles.info(2)-handles.Points_position(:,2);
%     
%     handles.Points_position(:,2) = handles.info(2)-handles.Points_position(:,2);
    Points_position = handles.Points_position;
    Points_color = handles.Points_color;

    save(handles.save_matfile, 'Points_position','Points_color');
end


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% redraw the image in the axes1(actually only add the filters and points according to
%% whether fliter=1  drawPoints=1 )
function [Newhandles]=redraw_image_inMainwin(filter,drawPoints,handles)
%global num;
global num_all;
 
if ( filter==1)
    %% if there are filters selected
    contents = cellstr(get(handles.popupmenu_filters_lowpass,'String')); 
    filter_name =  contents{get(handles.popupmenu_filters_lowpass,'Value')};

    if (strcmp(filter_name,'Gaussian'))     % Gaussian filter
        H = fspecial('gaussian', handles.hsize, handles.sigma);
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
    elseif (strcmp(filter_name,'Average'))       % average filter
        H = fspecial('average', handles.hsize);
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
    end
    contents = cellstr(get(handles.popupmenu_filters_highpass,'String')); 
    filter_name =  contents{get(handles.popupmenu_filters_highpass,'Value')};
    if (strcmp(filter_name,'Wiener'))
        handles.Grayim = wiener2(handles.Grayim,[handles.wsize handles.hsize]);
    elseif (strcmp(filter_name,'Sharp'))
        H = fspecial('unsharp');
        handles.Grayim = imfilter(handles.Grayim,H,'replicate');
    end
end
im = handles.Grayim;
scale = handles.scale;
set(handles.ih,'Cdata',im)
axes(handles.axes1);
hold on;
contents = cellstr(get(handles.popupmenu_color,'String'));
val_color=contents{get(handles.popupmenu_color,'Value')};   % get color
          
if (num_all >0) && (drawPoints==1)  % draw all the points
    for i=1: num_all   
        cx= handles.AllPoints_position(i,1); cy=handles.AllPoints_position(i,2);
        handles.cross(i) = plot(cx/scale,cy/scale,[handles.AllPoints_color(i) '*']);
    end

end
hold off;
set(handles.Points_num,'string', sprintf('Points Num: %d', num_all));  % show the Points' num in GUI
set(handles.Points_num, 'Value', num_all);

Newhandles = handles;





% --- Executes on button press in togglebutton_invertImage.
function togglebutton_invertImage_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_invertImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_invertImage
if handles.fileLoaded
 
%     button_state = get(hObject,'Value');
%     if button_state == get(hObject,'Max')
% 	% Toggle button is pressed, take appropriate action
%         a= size(handles.Grayim);  
%         white_mat = uint8(255.*ones(a));
%         im = white_mat-handles.Grayim;
%         set(handles.ih,'Cdata',im);
% %        set(hObject,'string','original image');
%     else
%         set(handles.ih,'Cdata',handles.Grayim);
% %        set(hObject,'string','invert image');
% 
%     end
%         
        im1 = get(handles.ih,'Cdata');
        
        a= size(im1);  
        white_mat = uint8(255.*ones(a));
        im = white_mat-im1;
        set(handles.ih,'Cdata',im);
        handles.Grayim =  im;
        
        
%%         Update handles structure
        guidata(hObject, handles);
        
    
end
      


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.fileLoaded
    MousePos = get (handles.axes1, 'CurrentPoint');
    set(handles.text_pos,'string',[ num2str(MousePos(1,1)) '    ' num2str(MousePos(1,2))]); %%p(1,1) is cols in image matrix(X coordinate)  p(1,2) is rows in image matrix(Y coordinate),matrix origin is on the left top

    api = iptgetapi(handles.hSP);  
    Vpart_pos= api.getVisibleLocation();

    pos = getpixelposition(handles.Main_im);
    scale = handles.scale;
    hightim = size(handles.originalim,1)/scale;
    widthim = size(handles.originalim,2)/scale;
     
    % judge whether the user clicks the image area
    if (MousePos(1,2) < (Vpart_pos(2)+min(pos(4),hightim))) && (MousePos(1,2) >Vpart_pos(2)) &&(MousePos(1,1) < (Vpart_pos(1)+min(pos(3),widthim))) && (MousePos(1,1) >Vpart_pos(1)) 
        setptr(gcf, 'crosshair');
    else
        setptr(gcf,'arrow');
    end
end

% % if handles.fileLoaded
% % obj_han=overobj('axes');
% % 
% % if obj_han
% %     set(gcf,'Pointer','crosshair');
% % else
% %     set(gcf,'Pointer','arrow');
% % end
% % end
% 
% obj_han=overobj('handles.Main_im');
% if ~isempty(obj_han) 
%     disp('hello'); 
% end;
