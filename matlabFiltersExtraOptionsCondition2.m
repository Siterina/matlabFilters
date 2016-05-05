function varargout = matlabFiltersExtraOptionsCondition2(varargin)
% MATLABFILTERSEXTRAOPTIONSCONDITION2 MATLAB code for matlabFiltersExtraOptionsCondition2.fig
%      MATLABFILTERSEXTRAOPTIONSCONDITION2, by itself, creates a new MATLABFILTERSEXTRAOPTIONSCONDITION2 or raises the existing
%      singleton*.
%
%      H = MATLABFILTERSEXTRAOPTIONSCONDITION2 returns the handle to a new MATLABFILTERSEXTRAOPTIONSCONDITION2 or the handle to
%      the existing singleton*.
%
%      MATLABFILTERSEXTRAOPTIONSCONDITION2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATLABFILTERSEXTRAOPTIONSCONDITION2.M with the given input arguments.
%
%      MATLABFILTERSEXTRAOPTIONSCONDITION2('Property','Value',...) creates a new MATLABFILTERSEXTRAOPTIONSCONDITION2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matlabFiltersExtraOptionsCondition2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matlabFiltersExtraOptionsCondition2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matlabFiltersExtraOptionsCondition2

% Last Modified by GUIDE v2.5 05-May-2016 20:09:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matlabFiltersExtraOptionsCondition2_OpeningFcn, ...
                   'gui_OutputFcn',  @matlabFiltersExtraOptionsCondition2_OutputFcn, ...
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


% --- Executes just before matlabFiltersExtraOptionsCondition2 is made visible.
function matlabFiltersExtraOptionsCondition2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matlabFiltersExtraOptionsCondition2 (see VARARGIN)

% Choose default command line output for matlabFiltersExtraOptionsCondition2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes1);
I = imread('condition2_equation.png');
imshow(I, [0, 2]);
% UIWAIT makes matlabFiltersExtraOptionsCondition2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matlabFiltersExtraOptionsCondition2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushButtonOk.
function pushButtonOk_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
omegacond2 = str2double(get(handles.editOmega, 'String'));
alphacond2 = str2double(get(handles.editAlpha, 'String'));
betacond2 = str2double(get(handles.editBeta, 'String'));
d11cond2 = str2double(get(handles.editD11, 'String'));
d22cond2 = str2double(get(handles.editD22, 'String'));
SoXcond2 = str2double(get(handles.editSoX, 'String'));
MoXcond2 = str2double(get(handles.editMoX, 'String'));
MoYcond2 = str2double(get(handles.editMoY, 'String'));
gamma = str2double(get(handles.editGamma, 'String'));

matlabFilters.initialConditions(0, -0.5, -1, 1.5, 0, 1, 0.5, 1, 0.1, 0.5, 0.6, ...
                                omegacond2, alphacond2, betacond2, d11cond2, d22cond2, SoXcond2, MoXcond2, MoYcond2, gamma);

close(handles.figure1);



function editOmega_Callback(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOmega as text
%        str2double(get(hObject,'String')) returns contents of editOmega as a double


% --- Executes during object creation, after setting all properties.
function editOmega_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to editAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAlpha as text
%        str2double(get(hObject,'String')) returns contents of editAlpha as a double


% --- Executes during object creation, after setting all properties.
function editAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBeta_Callback(hObject, eventdata, handles)
% hObject    handle to editBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBeta as text
%        str2double(get(hObject,'String')) returns contents of editBeta as a double


% --- Executes during object creation, after setting all properties.
function editBeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD11_Callback(hObject, eventdata, handles)
% hObject    handle to editD11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editD11 as text
%        str2double(get(hObject,'String')) returns contents of editD11 as a double


% --- Executes during object creation, after setting all properties.
function editD11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD22_Callback(hObject, eventdata, handles)
% hObject    handle to editD22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editD22 as text
%        str2double(get(hObject,'String')) returns contents of editD22 as a double


% --- Executes during object creation, after setting all properties.
function editD22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSoX_Callback(hObject, eventdata, handles)
% hObject    handle to editSoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSoX as text
%        str2double(get(hObject,'String')) returns contents of editSoX as a double


% --- Executes during object creation, after setting all properties.
function editSoX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMoX_Callback(hObject, eventdata, handles)
% hObject    handle to editMoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMoX as text
%        str2double(get(hObject,'String')) returns contents of editMoX as a double


% --- Executes during object creation, after setting all properties.
function editMoX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMoY_Callback(hObject, eventdata, handles)
% hObject    handle to editMoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMoY as text
%        str2double(get(hObject,'String')) returns contents of editMoY as a double


% --- Executes during object creation, after setting all properties.
function editMoY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGamma_Callback(hObject, eventdata, handles)
% hObject    handle to editGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGamma as text
%        str2double(get(hObject,'String')) returns contents of editGamma as a double


% --- Executes during object creation, after setting all properties.
function editGamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
