function varargout = water_filling(varargin)
% WATER_FILLING
%
% Visual representation of the water-filling (or water pouring) algorithm
% that allocates more (or less) bits and power to some subcarriers with
% larger (or smaller) SNR for % maximizing the channel capacity in OFDM systems.

% Giulio Marin
%
% giulio.marin@me.com
% 2013/03/20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @water_filling_OpeningFcn, ...
    'gui_OutputFcn',  @water_filling_OutputFcn, ...
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


% --- Executes just before water_filling is made visible.
function water_filling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to water_filling (see VARARGIN)


% Choose default command line output for water_filling
handles.output = hObject;

set(handles.Ptot_val,'String',round(get(handles.Ptot_slider,'Value')))
set(handles.SNR_val,'String',round(get(handles.SNR_slider,'Value')))

handles.Ptot=get(handles.Ptot_slider,'Value');
handles.SNR=get(handles.SNR_slider,'Value');
handles.ax1=handles.axes1;
handles.ax2=handles.axes2;
handles = update_plot(handles,'newChannel');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes water_filling wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = water_filling_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function SNR_slider_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SNR=get(hObject,'Value');
handles = update_plot(handles);
guidata(hObject, handles);
set(handles.SNR_val,'String',round(get(handles.SNR_slider,'Value')))

% --- Executes during object creation, after setting all properties.
function SNR_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Ptot_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Ptot_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.Ptot=get(hObject,'Value');
handles = update_plot(handles);
guidata(hObject, handles);
set(handles.Ptot_val,'String',round(get(handles.Ptot_slider,'Value')))


% --- Executes during object creation, after setting all properties.
function Ptot_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ptot_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Ptot_val_Callback(hObject, eventdata, handles)
% hObject    handle to Ptot_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ptot_val as text
%        str2double(get(hObject,'String')) returns contents of Ptot_val as a double


% --- Executes during object creation, after setting all properties.
function Ptot_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ptot_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_val_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR_val as text
%        str2double(get(hObject,'String')) returns contents of SNR_val as a double


% --- Executes during object creation, after setting all properties.
function SNR_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in channel_btn.
function channel_btn_Callback(hObject, eventdata, handles)
% hObject    handle to channel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = update_plot(handles,'newChannel');
guidata(hObject, handles);
