function varargout = PM_Analysis_AddParameter(varargin)
%PM_ANALYSIS_ADDPARAMETER M-file for PM_Analysis_AddParameter.fig
%      PM_ANALYSIS_ADDPARAMETER, by itself, creates a new PM_ANALYSIS_ADDPARAMETER or raises the existing
%      singleton*.
%
%      H = PM_ANALYSIS_ADDPARAMETER returns the handle to a new PM_ANALYSIS_ADDPARAMETER or the handle to
%      the existing singleton*.
%
%      PM_ANALYSIS_ADDPARAMETER('Property','Value',...) creates a new PM_ANALYSIS_ADDPARAMETER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to PM_Analysis_AddParameter_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PM_ANALYSIS_ADDPARAMETER('CALLBACK') and PM_ANALYSIS_ADDPARAMETER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PM_ANALYSIS_ADDPARAMETER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PM_Analysis_AddParameter

% Last Modified by GUIDE v2.5 19-Aug-2015 12:26:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PM_Analysis_AddParameter_OpeningFcn, ...
                   'gui_OutputFcn',  @PM_Analysis_AddParameter_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before PM_Analysis_AddParameter is made visible.
function PM_Analysis_AddParameter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for PM_Analysis_AddParameter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PM_Analysis_AddParameter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PM_Analysis_AddParameter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
    global xAddAnalysis
    %These global variables will be automatically updated in the main
    %script
    
    xAp = get(findobj('Tag','xanalysisParameter'),'Value');
    
    xT(1) = str2double(get(findobj('Tag','xetFrom'),'String'));
    xT(2) = str2double(get(findobj('Tag','xetTo'),'String'));
    
    xBp(1) = str2double(get(findobj('Tag','xetHigh'),'String'));
    xBp(2) = str2double(get(findobj('Tag','xetLow'),'String'));
    
    xGUI = get(findobj('Tag','xGUIselect'),'Value');
    
    xAddAnalysis = {xAp xT xBp xGUI};
    
    % close figure after done
    close gcbf;

function xetFrom_Callback(hObject, eventdata, handles)
% hObject    handle to xetFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetFrom as text
%        str2double(get(hObject,'String')) returns contents of xetFrom as a double


% --- Executes during object creation, after setting all properties.
function xetFrom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetTo_Callback(hObject, eventdata, handles)
% hObject    handle to xetTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetTo as text
%        str2double(get(hObject,'String')) returns contents of xetTo as a double


% --- Executes during object creation, after setting all properties.
function xetTo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetHigh_Callback(hObject, eventdata, handles)
% hObject    handle to xetHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetHigh as text
%        str2double(get(hObject,'String')) returns contents of xetHigh as a double


% --- Executes during object creation, after setting all properties.
function xetHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetLow2_Callback(hObject, eventdata, handles)
% hObject    handle to xetLow2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetLow2 as text
%        str2double(get(hObject,'String')) returns contents of xetLow2 as a double


% --- Executes during object creation, after setting all properties.
function xetLow2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetLow2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xanalysisParameter2.
function xanalysisParameter2_Callback(hObject, eventdata, handles)
% hObject    handle to xanalysisParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xanalysisParameter2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xanalysisParameter2


% --- Executes during object creation, after setting all properties.
function xanalysisParameter2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xanalysisParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in xGUIselect.
function xGUIselect_Callback(hObject, eventdata, handles)
    check = get(hObject,'Value');
    if check == 1; %(true)
        set(findobj('Tag','xetFrom'),'Enable','off');
        set(findobj('Tag','xetTo'),'Enable','off');
    else
        set(findobj('Tag','xetFrom'),'Enable','on');
        set(findobj('Tag','xetTo'),'Enable','on');
    end

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function xGUIselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xGUIselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in save2.
function save2_Callback(hObject, eventdata, handles)
    global xAddAnalysis
    %These global variables will be automatically updated in the main
    %script
    
    temp = get(findobj('Tag','listAnalysisOptions'),'String');
    
    for i = 1:size(temp,1)
        if i == 1
            temp2 = textscan(temp{1},'%f');
            temp2 = temp2{1}';
            xAddAnalysis = {temp2(1) [temp2(2) temp2(3)] [temp2(5) temp2(6)] temp2(4)};
        else
            temp2 = textscan(temp{i},'%f');
            temp2 = temp2{1}';
            xAddAnalysis(i,:) = {temp2(1) [temp2(2) temp2(3)] [temp2(5) temp2(6)] temp2(4)};
        end
    end
    % close figure after done
    close gcbf;


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on selection change in xanalysisParameter2.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to xanalysisParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xanalysisParameter2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xanalysisParameter2


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xanalysisParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetHigh2_Callback(hObject, eventdata, handles)
% hObject    handle to xetHigh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetHigh2 as text
%        str2double(get(hObject,'String')) returns contents of xetHigh2 as a double


% --- Executes during object creation, after setting all properties.
function xetHigh2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetHigh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetTo2_Callback(hObject, eventdata, handles)
% hObject    handle to xetTo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetTo2 as text
%        str2double(get(hObject,'String')) returns contents of xetTo2 as a double


% --- Executes during object creation, after setting all properties.
function xetTo2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetTo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xetFrom2_Callback(hObject, eventdata, handles)
% hObject    handle to xetFrom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xetFrom2 as text
%        str2double(get(hObject,'String')) returns contents of xetFrom2 as a double


% --- Executes during object creation, after setting all properties.
function xetFrom2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetFrom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listAnalysisOptions.
function listAnalysisOptions_Callback(hObject, eventdata, handles)
% hObject    handle to listAnalysisOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listAnalysisOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listAnalysisOptions


% --- Executes during object creation, after setting all properties.
function listAnalysisOptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listAnalysisOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddToList.
function AddToList_Callback(hObject, eventdata, handles)
    %These global variables will be automatically updated in the main
    %script
    
    xAp = get(findobj('Tag','xanalysisParameter2'),'Value');
    
    xT(1) = str2double(get(findobj('Tag','xetFrom2'),'String'));
    xT(2) = str2double(get(findobj('Tag','xetTo2'),'String'));
    
    xBp(1) = str2double(get(findobj('Tag','xetHigh2'),'String'));
    xBp(2) = str2double(get(findobj('Tag','xetLow2'),'String'));
    
    % Update list
    temp = get(findobj('Tag','listAnalysisOptions'),'String');
    if isempty(temp)
        temp = {num2str([xAp xT 0 xBp])};
        set(findobj('Tag','listAnalysisOptions'),'String',temp)
    else
        temp{end + 1,:} = num2str([xAp xT 0 xBp]);
        set(findobj('Tag','listAnalysisOptions'),'String',temp)
    end


% --- Executes on button press in RemoveFromList.
function RemoveFromList_Callback(hObject, eventdata, handles)

temp = get(findobj('Tag','listAnalysisOptions'),'String');
temp(get(findobj('Tag','listAnalysisOptions'),'Value'),:) = [];
set(findobj('Tag','listAnalysisOptions'),'String',temp);


% --- Executes during object creation, after setting all properties.
function xetLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xetLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function xanalysisParameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xanalysisParameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
