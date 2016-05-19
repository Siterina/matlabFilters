function varargout = matlabFiltersExtraOptionsCondition4(varargin)
% MATLABFILTERSEXTRAOPTIONSCONDITION4 MATLAB code for matlabFiltersExtraOptionsCondition4.fig
%      MATLABFILTERSEXTRAOPTIONSCONDITION4, by itself, creates a new MATLABFILTERSEXTRAOPTIONSCONDITION4 or raises the existing
%      singleton*.
%
%      H = MATLABFILTERSEXTRAOPTIONSCONDITION4 returns the handle to a new MATLABFILTERSEXTRAOPTIONSCONDITION4 or the handle to
%      the existing singleton*.
%
%      MATLABFILTERSEXTRAOPTIONSCONDITION4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATLABFILTERSEXTRAOPTIONSCONDITION4.M with the given input arguments.
%
%      MATLABFILTERSEXTRAOPTIONSCONDITION4('Property','Value',...) creates a new MATLABFILTERSEXTRAOPTIONSCONDITION4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matlabFiltersExtraOptionsCondition4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matlabFiltersExtraOptionsCondition4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matlabFiltersExtraOptionsCondition4

% Last Modified by GUIDE v2.5 20-May-2016 01:51:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matlabFiltersExtraOptionsCondition4_OpeningFcn, ...
                   'gui_OutputFcn',  @matlabFiltersExtraOptionsCondition4_OutputFcn, ...
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


% --- Executes just before matlabFiltersExtraOptionsCondition4 is made visible.
function matlabFiltersExtraOptionsCondition4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matlabFiltersExtraOptionsCondition4 (see VARARGIN)

% Choose default command line output for matlabFiltersExtraOptionsCondition4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes1);
I = imread('condition4_equation.png');
imshow(I);
truesize(handles.figure1);
% UIWAIT makes matlabFiltersExtraOptionsCondition4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matlabFiltersExtraOptionsCondition4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editA1_Callback(hObject, eventdata, handles)
% hObject    handle to editA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editA1 as text
%        str2double(get(hObject,'String')) returns contents of editA1 as a double


% --- Executes during object creation, after setting all properties.
function editA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editB1_Callback(hObject, eventdata, handles)
% hObject    handle to editB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editB1 as text
%        str2double(get(hObject,'String')) returns contents of editB1 as a double


% --- Executes during object creation, after setting all properties.
function editB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editC1_Callback(hObject, eventdata, handles)
% hObject    handle to editC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editC1 as text
%        str2double(get(hObject,'String')) returns contents of editC1 as a double


% --- Executes during object creation, after setting all properties.
function editC1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editC2_Callback(hObject, eventdata, handles)
% hObject    handle to editC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editC2 as text
%        str2double(get(hObject,'String')) returns contents of editC2 as a double


% --- Executes during object creation, after setting all properties.
function editC2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD1_Callback(hObject, eventdata, handles)
% hObject    handle to editD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editD1 as text
%        str2double(get(hObject,'String')) returns contents of editD1 as a double


% --- Executes during object creation, after setting all properties.
function editD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in pushButtonOk.
function pushButtonOk_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a1cond4 = str2double(get(handles.editA1, 'String'));
b1cond4 = str2double(get(handles.editB1, 'String'));
c1cond4 = str2double(get(handles.editC1, 'String'));
c2cond4 = str2double(get(handles.editC2, 'String'));
d1cond4 = str2double(get(handles.editD1, 'String'));
omegacond4 = str2double(get(handles.editOmega, 'String'));
MoXcond4 = str2double(get(handles.editMoX, 'String'));
SoXcond4 = str2double(get(handles.editSoX, 'String'));

matlabFilters.initialConditions(0, -0.5, -1, 1.5, 0, 1, 0.5, 1, 0.1, 0.5, 0.6, ...
                                0.1*3.1415, 2, 1, 0.1, 0.1, 0.58, 0.58, 10, -3, 1, ...
                                a1cond4, b1cond4, c1cond4, c2cond4, d1cond4, omegacond4, MoXcond4, SoXcond4);
close(handles.figure1);



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
