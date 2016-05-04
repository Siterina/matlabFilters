function varargout = matlabFiltersGUI(varargin)
% MATLABFILTERSGUI MATLAB code for matlabFiltersGUI.fig
%      MATLABFILTERSGUI, by itself, creates a new MATLABFILTERSGUI or raises the existing
%      singleton*.
%
%      H = MATLABFILTERSGUI returns the handle to a new MATLABFILTERSGUI or the handle to
%      the existing singleton*.
%
%      MATLABFILTERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATLABFILTERSGUI.M with the given input arguments.
%
%      MATLABFILTERSGUI('Property','Value',...) creates a new MATLABFILTERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matlabFiltersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matlabFiltersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matlabFiltersGUI

% Last Modified by GUIDE v2.5 05-May-2016 00:20:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matlabFiltersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @matlabFiltersGUI_OutputFcn, ...
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


% --- Executes just before matlabFiltersGUI is made visible.
function matlabFiltersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matlabFiltersGUI (see VARARGIN)

% Choose default command line output for matlabFiltersGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes matlabFiltersGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matlabFiltersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttonFilter.
function buttonFilter_Callback(hObject, eventdata, handles)
% hObject    handle to buttonFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
condition = 1;
buildTrajectory = false;
buildLAOF = true;
buildLFOS = true;
buildMx = false;

if (get(handles.radioButtonCondition1,'Value') == 1)
    condition = 1;
end
if (get(handles.radioButtonCondition2,'Value') == 1)
    condition = 2;
end
if (get(handles.radioButtonCondition3,'Value') == 1)
    condition = 3;
end
if(get(handles.checkBoxTrajectory, 'Value') == 1)
    buildTrajectory = true;
end
if(get(handles.checkBoxLAOF, 'Value') == 0)
    buildLAOF = false;
end
if(get(handles.checkBoxLFOS, 'Value') == 0)
    buildLFOS = false;
end
if(get(handles.checkBoxMx, 'Value') == 1)
    buildMx = true;
end

N = str2double(get(handles.editN, 'String'));
K = str2double(get(handles.editK, 'String'));
deltaTime = str2double(get(handles.editDeltaTime, 'String'));
trajectoryNumber = str2double(get(handles.editTrajectoryNumber, 'String'));

matlabFilters.main(condition, N, K, deltaTime, buildTrajectory, trajectoryNumber, buildLAOF, buildLFOS, buildMx);



function editN_Callback(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editN as text
%        str2double(get(hObject,'String')) returns contents of editN as a double


% --- Executes during object creation, after setting all properties.
function editN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDeltaTime_Callback(hObject, eventdata, handles)
% hObject    handle to editDeltaTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDeltaTime as text
%        str2double(get(hObject,'String')) returns contents of editDeltaTime as a double


% --- Executes during object creation, after setting all properties.
function editDeltaTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDeltaTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editK_Callback(hObject, eventdata, handles)
% hObject    handle to editK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editK as text
%        str2double(get(hObject,'String')) returns contents of editK as a double


% --- Executes during object creation, after setting all properties.
function editK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkBoxTrajectory.
function checkBoxTrajectory_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxTrajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxTrajectory



function editTrajectoryNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editTrajectoryNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrajectoryNumber as text
%        str2double(get(hObject,'String')) returns contents of editTrajectoryNumber as a double


% --- Executes during object creation, after setting all properties.
function editTrajectoryNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTrajectoryNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkBoxLAOF.
function checkBoxLAOF_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxLAOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxLAOF


% --- Executes on button press in checkBoxLFOS.
function checkBoxLFOS_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxLFOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxLFOS


% --- Executes on button press in checkBoxGAOF.
function checkBoxGAOF_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxGAOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxGAOF


% --- Executes on button press in checkBoxGFOS.
function checkBoxGFOS_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxGFOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxGFOS


% --- Executes on button press in checkBoxMx.
function checkBoxMx_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxMx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBoxMx
