function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates attenuazione new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to attenuazione new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates attenuazione new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 21-May-2016 01:04:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_OpeningFcn, ...
    'gui_OutputFcn',  @gui_OutputFcn, ...
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
end

% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.textNu,'String',round(get(handles.Nuslider,'Value')))
set(handles.textR,'String',[num2str(round(get(handles.Rslider,'Value'))), ' m'])

handles.ux=0;
handles.uy=0;
handles.R=get(handles.Rslider,'Value');
handles.Nu=get(handles.Nuslider,'Value');
handles.h1=handles.axes1;
handles.h2=handles.axes2;
handles.translation=get(handles.translider,'Value');
handles.B=3.84e6;
handles.SINRmin=1;
handles.varshadow=6;
handles.Pout = 0.1;
handles.a=3.8;
handles.Teff=2900;
handles.stampa=1;
handles.leg=[];
set(handles.editnote,'String','Modify the parameters');
set(handles.B_val,'String',handles.B)
set(handles.sinr_val,'String',handles.SINRmin)
set(handles.Pout_val,'String',handles.Pout)
set(handles.vars_val,'String',handles.varshadow)
set(handles.a_val,'String',handles.a)
set(handles.Teff_val,'String',handles.Teff)
cdma(handles);

% Update handles structure
guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on slider movement.
function translider_Callback(hObject, eventdata, handles)
% hObject    handle to translider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% data = get(handles.posizione,'UserData');
handles.translation=get(hObject,'Value');
handles = cdma(handles);
guidata(hObject, handles);
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function translider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to translider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have attenuazione light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
end

% --- Executes on slider movement.
function Rslider_Callback(hObject, eventdata, handles)
% hObject    handle to Rslider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.R=round(get(hObject,'Value'));
handles = cdma(handles);
guidata(hObject, handles);
set(handles.editnote,'String',handles.note);
set(handles.textR,'String',[num2str(round(get(handles.Rslider,'Value'))), ' m'])
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function Rslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rslider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have attenuazione light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on slider movement.
function Nuslider_Callback(hObject, eventdata, handles)
% hObject    handle to Nuslider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.Nu=round(get(handles.Nuslider,'Value'));
handles = cdma(handles);
guidata(hObject, handles);
set(handles.editnote,'String',handles.note);
set(handles.textNu,'String',round(get(handles.Nuslider,'Value')))
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function Nuslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nuslider (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have attenuazione light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


% --- Executes during object creation, after setting all properties.
function textR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textR (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end

% --- Executes during object creation, after setting all properties.
function textNu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textNu (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
end


% --- Executes on button press in posizione.
function posizione_Callback(hObject, eventdata, handles)
% hObject    handle to posizione (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2)
[handles.ux,handles.uy]=ginput(1);
handles = cdma(handles);
guidata(hObject, handles);
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end


% --- Executes on button press in creagrafico.
function creagrafico_Callback(hObject, eventdata, handles)
% hObject    handle to creagrafico (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.editnote,'String','Select a point to specify the direction');
axes(handles.axes2)
[handles.ux,handles.uy]=ginput(1);
handles = cdma(handles,'create');
guidata(hObject, handles);
set(handles.editnote,'String','Created plot Prmin(d/R)');
end



function sinr_val_Callback(hObject, eventdata, handles)
% hObject    handle to sinr_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sinr_val as text
%         returns contents of sinr_val as attenuazione double

handles.SINRmin=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end


% --- Executes during object creation, after setting all properties.
function sinr_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sinr_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have attenuazione white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Pout_val_Callback(hObject, eventdata, handles)
% hObject    handle to Pout_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pout_val as text
%        str2double(get(hObject,'String')) returns contents of Pout_val as attenuazione double
handles.Pout=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function Pout_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pout_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have attenuazione white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function a_val_Callback(hObject, eventdata, handles)
% hObject    handle to a_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_val as text
%        str2double(get(hObject,'String')) returns contents of a_val as attenuazione double
handles.a=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function a_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have attenuazione white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function B_val_Callback(hObject, eventdata, handles)
% hObject    handle to B_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_val as text
%        str2double(get(hObject,'String')) returns contents of B_val as attenuazione double
handles.B=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function B_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have attenuazione white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function vars_val_Callback(hObject, eventdata, handles)
% hObject    handle to vars_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vars_val as text
%        str2double(get(hObject,'String')) returns contents of vars_val as attenuazione double
handles.varshadow=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function vars_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vars_val (see GCBO)
% eventdata  reserved - to be defined in attenuazione future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have attenuazione white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Teff_val_Callback(hObject, eventdata, handles)
% hObject    handle to Teff_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Teff_val as text
%        str2double(get(hObject,'String')) returns contents of Teff_val as a double
handles.Teff=str2double(get(hObject,'String'));
handles = cdma(handles);
guidata(hObject, handles);
uicontrol(hObject)
set(handles.editnote,'String',handles.note);
set(handles.editprmin,'String',handles.Prmin);
end

% --- Executes during object creation, after setting all properties.
function Teff_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Teff_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


% Additional functions

function data = cdma(data, ~)

% TODO: impostare la soglia minima di PRmin nel grafico in base a Pn
% (potenza rumore). Prmin minimo è SINRmin*Pn

clc

if(nargin>1)
    d_R = linspace(0,1,30);
    ux = d_R*data.R*cos(atan2(data.uy,data.ux));
    uy = d_R*data.R*sin(atan2(data.uy,data.ux));
else
    d_R = 1;
    ux=data.ux;
    uy=data.uy;
end

data.note='';

Pd_R=[];
R=data.R;
Nu=round(data.Nu);
Pout=data.Pout;
SINRmin_dB=data.SINRmin;
sigma_s=data.varshadow;
sigma2_s=sigma_s^2;
B=data.B;
h1=data.h1;
h2=data.h2;
v=data.translation;
a=data.a;
Teff=data.Teff;
stampa=data.stampa;

%% PARAMETRI

K0 = (1/4/pi)^2;        % coefficiente di attenuazione
kB = 1.38e-23;          % costante di Boltzmann
G = 256;                % guadagno di processo
Prmin = 10.^((linspace(-152,-100,2e4)+v)/10);

alpha = 2/(a+2)/K0*Prmin*R^a;
beta = 2/(a+2)*Prmin*R^a/G;
k = log(10)/10;
SINRmin = 10^(SINRmin_dB/10);

%% Disegno celle

for j=1:length(d_R)
    
    axes(h2)
    plot(ux(j),uy(j),'bo','MarkerFaceColor','b')
    [~, ~, cx, cy] = hexagon(R,Nu);
        
    di = sqrt((cx-ux(j)).^2+(cy-uy(j)).^2);
    
    %% Shadowing singolo utente
    
    [m_s,var_s] = norm2log(0,k^2*sigma2_s);
    
    %% Potenza mean  trasmessa da una cella per utente
    
    mu_Pt = log(alpha);
    sigma2_Pt = k^2*sigma2_s;
    [m_Pt,var_Pt]=norm2log(mu_Pt,sigma2_Pt);
    
    %% Shadowing Nu utenti
    
    m_S = Nu*m_s;
    var_S = Nu*var_s;
    [mu_S,sigma2_S] = log2norm(m_S,var_S);
    
    %% Potenza mean  trasmessa da una cella
    
    mu_Pti = mu_S + log(alpha);
    sigma2_Pti = sigma2_S;
    [m_Pti,var_Pti] = log2norm(mu_Pti,sigma2_Pti);
    
    %% Potenza interferente ricevuta da ciascuna cella
    
    mu_PIi = k*10*log10(di'.^(-a)*beta)+mu_S;
    sigma2_PIi = sigma2_S+k^2*sigma2_s;
    [m_PIi,var_PIi] = norm2log(mu_PIi,sigma2_PIi);
    
    %% Potenza interferente totale
    
    mu_x = -k*a*10*log10(di)+mu_S;
    sigma2_x = sigma2_S+k^2*sigma2_s;
    
    [m_x,var_x] = norm2log(mu_x,sigma2_x);
    
    m_X = sum(m_x);
    var_X = sum(var_x);
    
    [mu_X, sigma2_X] = log2norm(m_X,var_X);
    
    mu_PI = mu_X+log(beta);
    sigma2_PI = sigma2_X;
    [m_PI,var_PI] = norm2log(mu_PI,sigma2_PI);
    
    %% Rumore termico
    
    Pn = kB*Teff*B/G;
    
    %% Grafico
    
    x = Prmin./SINRmin-Pn;
    
    x(x<0) = 0; % valori negativi hanno probabilità 0
    p = 1-normcdf(log(x),mu_PI,sigma2_PI);
    x1 = log(x);
    x2 = norminv(1-Pout,mu_PI,sigma2_PI);
    
    axes(h1);
    plot(10*log10(Prmin),p,10*log10([min(Prmin),max(Prmin)]),[Pout,Pout],'r--','linewidth',2)
    legend('current value','minimum value')
    xlabel('P_{Rmin} dB')
    ylabel('P_{out}')
    title('Outage probability')
    axis([10*log10([min(Prmin),max(Prmin)]) 0 1])
    grid on
    
    %% Risultati

    i = floor(length(Prmin)/2);
    if ((min(x1)<min(x2) && max(x1)>max(x2)) || (min(x2)<min(x1) && max(x2)>max(x1)))
        %     sprintf('SI\n')];
        i=find(abs(p-Pout)==min(abs(p-Pout)),1);
        Pd_R=[Pd_R, Prmin(i)];
        data.Prmin = Prmin(i);
        %     sprintf('(errore probabilità ');
        if (abs(Pout-p(i))<0.05)
            %         sprintf('inferiore al 5 percento)\n');
        else
            %         sprintf('superiore al 5 percento)\n');
            %         sprintf('\n!!! Restringi il campo di ricerca della potenza !!!\n\n');
        end
    else
        data.Prmin=0;
        if(min(p)<Pout)
            Pd_R=[Pd_R, 0];
            data.note='Prmin smaller than current limit';
        else
            Pd_R=[Pd_R, inf];
            if(Nu==0)
                data.note='Interferenza nulla';
                data.Prmin=SINRmin*Pn;
            else
                data.note='Prmin higher than current limit';
            end
        end
        
    end
    
end

if(nargin>1)
    data.leg{length(data.leg)+1}=['Nu=',num2str(Nu)];
    figure(1)
    set(1,'DefaultAxesLineStyleOrder',{'--o'})
    hold all
    plot(d_R,10*log10(Pd_R))
    legend(data.leg)
    xlabel('d/R')
    ylabel('P_{Rmin} dB')
    title('Minimum received power as a function of the distance')
    hold off
end

%% Stampa i valori

if(stampa)
    
    fprintf('\n-------------------------------------------\n')
    fprintf('             PARAMETERS                    ')
    fprintf('\n-------------------------------------------\n')
    
    
    fprintf('Number of users                %d\n',Nu)
    fprintf('Process gain                   %d\n',G)
    fprintf('SINRmin                        %d dB\n',SINRmin_dB)
    fprintf('Outage probability             %g\n',Pout)
    fprintf('Attenuation coefficient        %g\n',a)
    fprintf('Standard deviation shadowing   %d dB\n',sigma_s)
    fprintf('Band                           %g Hz\n',B)
    fprintf('Effective temperature          %d K\n',Teff)
    
    fprintf('\n-------------------------------------------\n')
    fprintf('              RESULTS                 ')
    fprintf('\n-------------------------------------------\n')
    
    
    fprintf('\n* Shadowing single user\n\n');
    fprintf('dB\n')
    fprintf('\t mean  \t\t %.3g\n',0);
    fprintf('\t variance\t %.3g\n\n',sigma2_s)
    
    fprintf('linear\n')
    fprintf('\t mean  \t\t %.3f\n',m_s);
    fprintf('\t variance\t %.3g\n',var_s)
    
    fprintf('\n* Shadowing Nu users\n\n');
    fprintf('dB\n')
    fprintf('\t mean  \t\t %.3g\n',mu_S/k);
    fprintf('\t variance\t %.3g\n\n',sigma2_S/k^2)
    
    fprintf('linear\n')
    fprintf('\t mean  \t\t %.3g\n',m_S);
    fprintf('\t variance\t %.3g\n',var_S)
    
    fprintf('\n* Transmitted poer per user\n\n');
    fprintf('dB\n')
    fprintf('\t mean  \t\t %.3g\n',mu_Pt(i)/k);
    fprintf('\t variance\t %.3g\n\n',sigma2_Pt/k^2)
    
    fprintf('linear\n')
    fprintf('\t mean  \t\t %.3g\n',m_Pt(i));
    fprintf('\t variance\t %.3g\n',var_Pt(i))
    
    
    fprintf('\n* Transmitted poer per  cell\n\n');
    fprintf('dB\n')
    fprintf('\t mean  \t\t %.3g\n',mu_Pti(i)/k);
    fprintf('\t variance\t %.3g\n\n',sigma2_Pti/k^2);
    
    fprintf('linear\n')
    fprintf('\t mean  \t\t %.3g\n',m_Pti(i));
    fprintf('\t variance\t %.3e\n',var_Pti(i))
    
    
    fprintf('\n* Interferent power of the cell i\n\n');
    fprintf('dB\n')
    fprintf('\t mean  %d \t %.3g\n',[1:6;(mu_PIi(:,i)/k)']);
    fprintf('\t variance\t %.3g\n\n',sigma2_PIi/k^2)
    
    fprintf('linear\n')
    fprintf('\t mean  %d \t %.3g\n',[1:6;m_PIi(:,i)']);
    fprintf('\t variance %d\t %.3g\n',[1:6;var_PIi(:,i)']);
    
    fprintf('\n* Noise power\n\n');
    fprintf('dB\n')
    fprintf('\t %.3g\n',10*log10(Pn));
    
    fprintf('linear\n')
    fprintf('\t %.3g\n',Pn);
    
    fprintf('\n*Noise power + total interference\n\n');
    fprintf('dB\n')
    fprintf('\t mean  \t\t %.3g\n',mu_PI(i)/k);
    fprintf('\t variance\t %.3g\n\n',sigma2_PI/k^2)
    
    fprintf('linear\n')
    fprintf('\t mean  \t\t %.3g\n',m_PI(i));
    fprintf('\t variance\t %.3g\n',var_PI(i));
    
    fprintf('\n* Minimum received power\n\n');
    fprintf('dB\n')
    fprintf('\t%.0f\n\n',10*log10(data.Prmin))
    
    fprintf('linear\n')
    fprintf('\t%.4g\n\n',data.Prmin)
    
end
end

function [x0,y0,cx,cy] = hexagon(r, Nu)
%HEXAGON Hexagonal cell structure.
%   HEXAGON(r) generates a 7-reuse hexagonal cell structure with radius r.
%
%   [x0,y0] 
%   [cx,cy] column vectors of the centers of the hexagons, starting from
%   the top one and proceeding clockwise
 
    hold on
    t=linspace(0,2*pi,7);               % Central hexagon
    x0=r*cos(t);
    y0=r*sin(t);
    plot(x0,y0);
    plot(0,0,'ro');
    grid on
    axis equal

    x1=x0;                              % Upper hex
    y1=y0+r*sqrt(3);
    plot(x1,y1);
    plot(0,r*sqrt(3),'ko');

    x2=x0;                              % Lower hex
    y2=y0-r*sqrt(3);
    plot(x2,y2);
    plot(0,-r*sqrt(3),'ko');

    x3=x0+3*r/2;                        % Right lower hex
    y3=y0-sqrt(3)*r/2;
    plot(x3,y3);
    plot(3*r/2,-sqrt(3)*r/2,'ko');

    x4=x0+3*r/2;                        % Right upper hex
    y4=y0+sqrt(3)*r/2;
    plot(x4,y4);
    plot(3*r/2,sqrt(3)*r/2,'ko');

    x5=x0-3*r/2;                        % Left upper hex
    y5=y0+sqrt(3)*r/2;
    plot(x5,y5);
    plot(-3*r/2,sqrt(3)*r/2,'ko');

    x6=x0-3*r/2;                     % Left lower hex
    y6=y0-sqrt(3)*r/2;
    plot(x6,y6);
    plot(-3*r/2,-sqrt(3)*r/2,'ko'); 
    
    cx = [0,3*r/2,3*r/2,0,-3*r/2,-3*r/2];
    cy = [sqrt(3)*r,sqrt(3)*r/2,-sqrt(3)*r/2,-sqrt(3)*r,-sqrt(3)*r/2,sqrt(3)*r/2];

    % Draw users
    if false
        for i=1:6
            nx=-r+2*r*rand(1,Nu);
            ny=-r+2*r*rand(1,Nu);
            plot(nx+cx(i),ny+cy(i),'kx')
        end
    end
    axis([-3*r 3*r -3*r 3*r])
    hold off
    
end

function [m_v,var_v] = log2norm(m_s,var_s)
%LOG2NORM Returns mean and variance of the associated normal
%random variable.
%
%   s=exp(v) and v~N(m_v,var_v)

var_v = log(1+var_s./(m_s.^2));
m_v = log(m_s)-var_v/2;
end

function [m_s,var_s] = norm2log(m_v,var_v)
%NORM2LOG Returns mean and variance of the associated lognormal
%random variable.
%
%   s=exp(v) and v~N(m_v,var_v)

m_s = exp((2*m_v+var_v)/2);
var_s = m_s.^2.*(exp(var_v)-1);
end
