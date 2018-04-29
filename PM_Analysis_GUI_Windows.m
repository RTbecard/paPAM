function varargout = PM_Analysis_GUI_Windows(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PM_Analysis_GUI_Windows_OpeningFcn, ...
                   'gui_OutputFcn',  @PM_Analysis_GUI_Windows_OutputFcn, ...
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


% --- Executes just before PM_Analysis_GUI_Windows is made visible.
function PM_Analysis_GUI_Windows_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = PM_Analysis_GUI_Windows_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

 %% List saved profiles
    profile{1} = 'default';
    
    %% Load selected profile
    
    pathCalib = '';
    calibFile = '';
    pathFile = '';
    pathTone = '';
    fileName = '';
    output = '';
    
    try
        if ispc()
           load(['Profiles\' profile{1}]); 
        else
           load(['Profiles/' profile{1}]);
        end
        
    catch
       return; 
    end
    
    try %#ok<TRYNC>
        %% Update UI with loaded variables
        set(findobj('tag','Threshold'),'String',num2str(Threshold));
        set(findobj('tag','FileType'),'Value',FileType);
        set(findobj('tag','analysisType'),'Value',AnalysisType);
        set(findobj('tag','analysisParameter'),'Value',AnalysisParameter);
        set(findobj('tag','SampleRate'),'String',num2str(SampleRate));

        set(findobj('tag','etTo'),'String',num2str(TimeStart));
        set(findobj('tag','etFrom'),'String',num2str(TimeEnd));
        set(findobj('tag','etHigh'),'String',num2str(BandpassHigh));
        set(findobj('tag','etLow'),'String',num2str(BandPassLow));
        set(findobj('tag','PSDHz'),'String',num2str(PSDHz));

        set(findobj('tag','customTimeStamp'),'Value',timestampOption);

        set(findobj('tag','waveform'),'Value',WaveForm);
        set(findobj('tag','spectorgram'),'Value',Spectogram);
        set(findobj('tag','writeToFile'),'Value',WriteToFile);
        set(findobj('tag','preserveMemory'),'Value',preserveMemory);
        set(findobj('tag','median'),'Value',medianX);
        set(findobj('tag','saveWindowData'),'Value',setWindowData);

        set(findobj('tag','PVL'),'Value',PVL);
        set(findobj('tag','PAL'),'Value',PAL);

        set(findobj('tag','GUIselect'),'Value',GUISelect);

        set(findobj('tag','GainA'),'String',num2str(GainA));
        set(findobj('tag','GainB'),'String',num2str(GainB));
        set(findobj('tag','GainC'),'String',num2str(GainC));
        set(findobj('tag','GainD'),'String',num2str(GainD));
        set(findobj('tag','GainRef'),'String',num2str(GainR));

        set(findobj('tag','windowType'),'Value',WindowType);
        set(findobj('tag','windowLength'),'String',num2str(WindowLength));
        set(findobj('tag','overlap'),'String',num2str(Overlap));

        set(findobj('tag','figureSizeX'),'String',num2str(FigureSize(1)));
        set(findobj('tag','figureSizeY'),'String',num2str(FigureSize(2)));
    end
    
    fileName = '';
    pathFile = '';
    
    try
       save('pathCalib.mat','pathCalib','calibFile','pathFile','pathTone','fileName','output');
    catch
    end
    
    updateFilePaths;

function analysisType_Callback(hObject, eventdata, handles)
    if get(findobj('Tag','analysisType'),'Value') == 1;
        set(findobj('Tag','IgnoreFirstSeconds'),'Enable','off');
        set(findobj('Tag','Threshold'),'Enable','off');
        set(findobj('Tag','ThresholdWait'),'Enable','off');
        set(findobj('Tag','SelectAnalysisFile'),'String','Select File');
        set(findobj('Tag','writeToFile'),'Enable','on');
    else
        set(findobj('Tag','IgnoreFirstSeconds'),'Enable','on');
        set(findobj('Tag','Threshold'),'Enable','on');
        set(findobj('Tag','ThresholdWait'),'Enable','on');
        set(findobj('Tag','SelectAnalysisFile'),'String','Select Folder');

        set(findobj('Tag','etTo'),'Enable','on');
        set(findobj('Tag','etFrom'),'Enable','on');

        set(findobj('Tag','writeToFile'),'Enable','off');
        set(findobj('Tag','writeToFile'),'Value',1);
    end

    updateFilePaths;


function buttonPath_Callback(hObject, eventdata, handles)
    path = uigetdir;
    tPath = findobj('tag','tPath');
    set(tPath,'String',path);

function selectAnalysisFile_Callback(hObject, eventdata, handles)
    
    pathCalib = '';
    pathFile = '';
    fileName = '';
    pathTone = '';
    calibFile = '';
    outputfolder = '';
    
    try
       load('pathCalib.mat');
    catch
    end

    if get(findobj('Tag','analysisType'),'Value') == 1;
        [fileName, pathFile] = uigetfile({'*.wav' '.wav';'*.csv' '.csv';'*.mp3' '.mp3'},'Select File',[pathFile fileName]); %#ok<ASGLU>
    else
        pathFile = uigetdir(pathFile,'Select Folder');
        fileName = '';
    end
    
    if pathFile > 0
        save('pathCalib.mat','pathCalib','pathFile','pathTone','fileName','calibFile','outputfolder');
        updateFilePaths;
    end
    


function FileType_Callback(hObject, eventdata, handles)
    if get(findobj('Tag','FileType'),'Value') < 3;
        set(findobj('Tag','SampleRate'),'Enable','on');
        set(findobj('Tag','SampleRate'),'BackgroundColor',[0.831 0.816 0.784]);

    else
        set(findobj('Tag','SampleRate'),'Enable','off');
        set(findobj('Tag','SampleRate'),'BackgroundColor',[0.831 0.816 0.784]);
    end

function selectcalibrationFile_Callback(hObject, eventdata, handles)

   pathCalib = '';
   pathFile = '';
   pathTone = '';
   fileName = '';
   calibFile = '';
   outputfolder = '';
    try
        load('pathCalib.mat');
    catch
    end
    
    if length(pathCalib < 2) %#ok<ISMT>
        pathCalib = '';
    end
    
    [calibFile, pathCalib] = uigetfile({'*.csv' '.csv comma separated file'},pathCalib); %#ok<ASGLU>
    
    save('pathCalib.mat','pathCalib','pathFile','pathTone','fileName','calibFile','outputfolder');
    updateFilePaths;

% --- Executes on button press in Proccess.
function Proccess_Callback(hObject, eventdata, handles)
    Threshold = str2double(get(findobj('tag','Threshold'),'String'));
    ThresholdWait = str2double(get(findobj('tag','ThresholdWait'),'String'));
    AutoDetectDelay = str2double(get(findobj('tag','IgnoreFirstSeconds'),'String'));
    FileType = get(findobj('tag','FileType'),'Value');
    AnalysisType = get(findobj('tag','analysisType'),'Value');
    AnalysisParameter = get(findobj('tag','analysisParameter'),'Value');
    
    SampleRate = str2double(get(findobj('tag','SampleRate'),'String'));
    
    TimeStart = str2double(get(findobj('tag','etTo'),'String'));
    TimeEnd = str2double(get(findobj('tag','etFrom'),'String'));
    BandpassHigh = str2double(get(findobj('tag','etHigh'),'String'));
    BandPassLow = str2double(get(findobj('tag','etLow'),'String'));

    saveWindowData = get(findobj('tag','saveWindowData'),'Value');
    
    timestampOption = get(findobj('tag','timeStamp'),'Value');
    timestamp = get(findobj('tag','customTimeStamp'),'String');
    
    %% Analysis options
    WaveForm = get(findobj('tag','waveform'),'Value');
    Spectogram = get(findobj('tag','spectorgram'),'Value');
    PVL = get(findobj('tag','PVL'),'Value');
    PAL = get(findobj('tag','PAL'),'Value');
    percentiles = get(findobj('tag','median'),'Value');
    
    WriteToFile = get(findobj('tag','writeToFile'),'Value');
    preserveMemory = get(findobj('tag','preserveMemory'),'Value');
    
    GUISelect = get(findobj('tag','GUIselect'),'Value');
    
    % Consistancy analysis parameters
    consistencyP = str2double(get(findobj('Tag','consistencyP'),'String'));
    consistencyV = str2double(get(findobj('Tag','consistencyV'),'String'));
    consistencyA = str2double(get(findobj('Tag','consistencyV'),'String'));
    
    %% Recorder gain
    GainA = str2double(get(findobj('tag','GainA'),'String'));
    GainB = str2double(get(findobj('tag','GainB'),'String'));
    GainC = str2double(get(findobj('tag','GainC'),'String'));
    GainD = str2double(get(findobj('tag','GainD'),'String'));
    
    WindowType = get(findobj('tag','windowType'),'Value');
    temp = get(findobj('tag','windowType'),'String');
    WindowType = temp{WindowType};
    
    singleChan = get(findobj('tag','singleChan'),'Value');
    
    if get(findobj('tag','sampleratewindow'),'Value') == 0
        WindowLength = str2double(get(findobj('tag','windowLength'),'String'));
    else
        WindowLength = 0;
    end
    Overlap = str2double(get(findobj('tag','overlap'),'String'));
    
    %% Figure options
    FigureSize(1) = str2double(get(findobj('tag','figureSizeX'),'String'));
    FigureSize(2) = str2double(get(findobj('tag','figureSizeY'),'String'));
    fontSize = str2double(get(findobj('tag','fontSize'),'String'));
    fontType = get(findobj('tag','fontType'),'Value');
    fonts = listfonts;
    fontType = fonts{fontType};
    
    publishableFigures = get(findobj('tag','publishableFigures'),'Value');

    colormapVal = get(findobj('tag','colormapH'),'Value');
    colormaps = get(findobj('tag','colormapH'),'String');
    colormapVal = colormaps{colormapVal};

    autoSelectColorbar = get(findobj('Tag','autoSelectColorbar'),'Value');
    colorLp = str2double(get(findobj('Tag','colorLp'),'String'));
    colorHp = str2double(get(findobj('Tag','colorHp'),'String'));
    colorLv = str2double(get(findobj('Tag','colorLv'),'String'));
    colorHv = str2double(get(findobj('Tag','colorHv'),'String'));
    colorLa = str2double(get(findobj('Tag','colorLa'),'String'));
    colorHa = str2double(get(findobj('Tag','colorHa'),'String'));
    
    AutoSelectFreqRange = get(findobj('Tag','AutoSelectFreqRange'),'Value');
    freqMin = str2double(get(findobj('Tag','freqMin'),'String'));
    freqMax = str2double(get(findobj('Tag','freqMax'),'String'));
    
    calibUnitPressure = get(findobj('Tag','CalibUnitPressure'),'Value');
    calibUnitVelocity = get(findobj('Tag','CalibUnitVelocity'),'Value');
    
    %% Load Paths
    pathCalib = '';
    pathFile = '';
    pathTone = '';
    fileName = '';
    calibFile = '';
    outputfolder = '';

    try
        load('pathCalib.mat');
    end
    
    %% Check for errors in user input
    AnalysisType = get(findobj('tag','analysisType'),'Value');
    if ispc
        temp = [pathCalib '\'];
        temp2 = [pathFile '\'];
    else
        temp = [pathCalib '/'];
        temp2 = [pathFile '/'];
    end
    run = 1;
    if exist([temp calibFile]) == 0 || isempty(calibFile)
       msgbox('No device calibration file specified, or file does not exist.','Select Calibration File!');
       return;
    end
    if AnalysisType == 1
       filecheck = exist([temp2 fileName]);
       if isempty(fileName) == 1
        filecheck = 0;
       end
    else
       filecheck = exist(temp2); %#ok<*EXIST>
    end
    
    if filecheck == 0
        msgbox('No input audio file specified, or file does not exist.','Select Input');
        return;
    end
    if exist(outputfolder) == 0
       msgbox('Invalid output folder selected.  Results will be saved to the current working directory.','Select Output Folder','modal');
       uiwait;
       outputfolder = pwd;
    end
    if (PVL + PAL == 0) && singleChan == 1
        beep;
        msgbox('You must select either velocity or acceleration as output.','Select Output','modal'); 
        return;
    end

    %% Execute analysis
    PM_Analysis_DataCrawler({pathFile fileName},{pathCalib calibFile},AnalysisType,AnalysisParameter,Threshold,ThresholdWait,...
        GUISelect,[TimeStart TimeEnd],[BandpassHigh BandPassLow],{WaveForm Spectogram WriteToFile PVL PAL outputfolder...
        singleChan percentiles saveWindowData preserveMemory},...
        [GainA GainB GainC GainD],SampleRate,{timestampOption timestamp},{WindowType WindowLength Overlap},...
        {FigureSize,fontSize,fontType,publishableFigures,colormapVal,autoSelectColorbar,[colorLp colorHp;colorLv colorHv;colorLa colorHa]},...
        [consistencyP consistencyV consistencyA],{AutoSelectFreqRange freqMin freqMax}, {calibUnitPressure calibUnitVelocity});
    disp('Analysis complete!');
    beep;            


% --- Executes on button press in selectCalibrationRecording.
function selectCalibrationRecording_Callback(hObject, eventdata, handles)
    % Load previously saved paths
    pathCalib = ''; %#ok<*NASGU>
    pathFile = '';
    pathTone = '';
    fileName = '';
    calibFile = '';
    outputfolder = '';
    try
        load('pathCalib.mat');
    catch
    end
    
    answer = inputdlg({'Reference Voltage (Voltage peak-to-peak)'},'Channel Calibration',1,{'0'});
    
    %Get calibration voltage 
    refVoltage = str2double(answer);
    
    if refVoltage == 0;
        msgbox('Reference voltage must be greater than 0.');
        return;
    end
    
    % Select device channel
    selection = listdlg('PromptString','Select recorder channel to calibrate','SelectionMode','single','ListString',{'X' 'Y' 'Z' 'Hydrophone'},'ListSize',[200 150]);
    %Select calibration tone
    [file, pathTone] = uigetfile({'*.wav','wav file'},'Select Calibration Recording',pathTone);
    % if user cancels file select, exit function
    if length(pathTone) < 2;
        return
    end
        
    save('pathCalib.mat','pathCalib','pathFile','pathTone','fileName','calibFile','outputfolder');
    channels = {'GainA' 'GainB' 'GainC' 'GainD'};
   
    % Load calibration tone into memory
    if ispc
        if  verLessThan('matlab', '8.1');
            [y,fs] = wavread([pathTone '\' file],'double'); 
        else
            [y,fs] = audioread([pathTone '\' file],'double'); 
        end
    else
        if  verLessThan('matlab', '8.1');
            [y,fs] = wavread([pathTone '/' file],'double'); 
        else
            [y,fs] = audioread([pathTone '/' file],'double'); 
        end
    end
  
    % if multiple channels in wav, ask user to select correct channel
    if size(y,2) == 1;
        channelwavS = 1;
    else    
        for i=1:size(y,2)
            channelWav{i} = num2str(i);
        end
        channelwavS = listdlg('PromptString','Select channel of calibration recording','SelectionMode','single','ListString',channelWav,'ListSize',[200 150]);
    end
    
    %Extract selected channel from recording
    y = y(:,channelwavS);
    
    figure(1);
    %plot amplitude vs seconds
    plot((1:(size(y,1)))./fs,y);
    grid on;
    xlabel('Seconds');
    ylabel('Amplitude');
    msg = msgbox('Select the start of the calibration tone');
    uiwait(msg);
    figure(1);
    [Sx(1),~] = ginput(1);
    maxAmp = max(abs(y));
    hold on
    plot([Sx(1) Sx(1)],[-maxAmp maxAmp],':r','LineWidth',2);
    hold off
    msg = msgbox('Select the end of the calibration tone');
    uiwait(msg);
    figure(1);
    [Sx(2),~] = ginput(1);
    close 1;
    
    disp(['start sample:' num2str(Sx(1))])
    disp(['end sample:'  num2str(Sx(2))])
    
    %Return trimmed calibration data
    trim = y(uint32(Sx(1)*fs):uint32(Sx(2)*fs));
    
    %get peak to peak value from recording
    calibPeakToPeak = max(trim)-min(trim);
    
    set(findobj('tag',channels{selection}),'String',num2str(refVoltage/calibPeakToPeak));

    

% --- Executes on button press in GUIselect.
function GUIselect_Callback(hObject, eventdata, handles)
    if get(hObject,'Value') == 1;
        set(findobj('Tag','etFrom'),'Enable','off');
        set(findobj('Tag','etTo'),'Enable','off');
    else
        set(findobj('Tag','etFrom'),'Enable','on');
        set(findobj('Tag','etTo'),'Enable','on');
    end

function logo_CreateFcn(hObject, eventdata, handles)
    [x,~,A] = imread('logo.png');
    h = imshow(x);
    axis image;
    set(h, 'AlphaData', A);

function updateFilePaths
    pathCalib = '';
    pathFile = '';
    pathTone = '';
    fileName = '';
    calibFile = '';
    outputfolder = '';
    try
        load('pathCalib.mat');
    catch
    end
    
    set(findobj('tag','DeviceCalibrationFile'),'String',calibFile);

    if length(outputfolder) > 1
        if ispc
            temp = textscan(outputfolder,'%s','delimiter','\\');
        else
            temp = textscan(outputfolder,'%s','delimiter','/');
        end
        temp = temp{1};
        temp = temp{end};
    else
        temp = '';
    end
    
    set(findobj('tag','outputfolder'),'String',temp);
    
    tPath = findobj('tag','tPath');

    AnalysisType = get(findobj('tag','analysisType'),'Value');
    if AnalysisType == 1
        set(tPath,'String',fileName);
    else
        if length(pathFile) > 1
            if ispc
                temp = textscan(pathFile,'%s','delimiter','\\');
            else
                temp = textscan(pathFile,'%s','delimiter','/');
            end
            temp = temp{1};
            temp = temp{end};
        else
            temp = '';
        end
        set(tPath,'String',temp);
    end


% --- Executes on button press in sampleratewindow.
function sampleratewindow_Callback(hObject, eventdata, handles)
    
    value = get(hObject,'Value');
    if value == 1;
        set(findobj('Tag','windowLength'),'Enable','off');
    else
        set(findobj('Tag','windowLength'),'Enable','on');
    end
    

% --- Executes on button press in setDefault.
function setDefault_Callback(hObject, eventdata, handles)
    %% Enter profile name 
    saveNewProfile('default');
    
    msgbox('The current settings have been saved as the default','Settings');

function Save_Callback(hObject, eventdata, handles)

    %% Enter profile name 
    profileName = inputdlg({'Enter name of profile:'},'Save profile',1,{'profile name'});
    
    saveNewProfile(profileName{1});

function Load_Callback(hObject, eventdata, handles)
    %% List saved profiles
    warning('off','all')
    mkdir('Profiles');
    warning('on','all')
    profileNames = dir('Profiles');
    for i=1:length(profileNames)
       profiles{i} = profileNames(i).name; 
    end

    selection = listdlg('ListString',profiles,'SelectionMode','single');
    profile = profiles(selection);
    
    if isempty(selection)
        return
    end
    
    %% Load file paths
    pathFile = '';
    pathTone = '';
    fileName = '';
    
    %% Load selected profile
    if ispc()
       ProfileLocation = (['Profiles\' profile{1}]); 
    else
       ProfileLocation = (['Profiles/' profile{1}]);
    end
    
    loadProfile(ProfileLocation);
    
    try
       save('pathCalib.mat','pathCalib','calibFile','pathFile','pathTone','fileName','outputfolder');
    end
    updateFilePaths;
   
function selectoutputfolder_Callback(hObject, eventdata, handles)
	pathCalib = '';
    pathFile = '';
    pathTone = '';
    fileName = '';
    calibFile = '';
    outputfolder = '';
    try
        load('pathCalib.mat');
    catch
    end
    
    if length(outputfolder) < 2
        outputfolder = pwd;
    end
    
    try
        [outputfolder] = uigetdir(outputfolder);
    catch
    end
    
    if length(outputfolder) < 2
        outputfolder = pwd;
    end
    
    save('pathCalib.mat','pathCalib','pathFile','pathTone','fileName','calibFile','outputfolder');
    updateFilePaths;


function pushbutton36_Callback(hObject, eventdata, handles)

if ispc();
    try
        winopen('Profiles');
    catch
        disp('No profiles folder');
    end
else
    disp('This button only works on windows platforms.  For other operating systems, you can find your saved profiles in the same directory where paPAM is saved.');
end

function fontType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
fonts = listfonts;
idx = find(strcmp('Helvetica',fonts));

if isempty(idx)
    idx = 1;
end

set(hObject,'String',fonts);
set(hObject,'Value',idx);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in autoSelectColorbar.
function autoSelectColorbar_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
       set(findobj('Tag','colorLp'),'Enable','off');
       set(findobj('Tag','colorHp'),'Enable','off');
       set(findobj('Tag','colorLv'),'Enable','off');
       set(findobj('Tag','colorHv'),'Enable','off');
       set(findobj('Tag','colorLa'),'Enable','off');
       set(findobj('Tag','colorHa'),'Enable','off');
    else
       set(findobj('Tag','colorLp'),'Enable','on');
       set(findobj('Tag','colorHp'),'Enable','on');
       set(findobj('Tag','colorLv'),'Enable','on');
       set(findobj('Tag','colorHv'),'Enable','on');
       set(findobj('Tag','colorLa'),'Enable','on');
       set(findobj('Tag','colorHa'),'Enable','on');
    end
    
function loadProfile(ProfileLocation)
%    try
        load(ProfileLocation);

        %% Update UI with loaded variables
        set(findobj('tag','Threshold'),'String',num2str(Threshold));
        set(findobj('tag','FileType'),'Value',FileType);
        set(findobj('tag','analysisType'),'Value',AnalysisType);
        set(findobj('tag','analysisParameter'),'Value',AnalysisParameter);
        set(findobj('tag','SampleRate'),'String',num2str(SampleRate));

        set(findobj('tag','etTo'),'String',num2str(TimeStart));
        set(findobj('tag','etFrom'),'String',num2str(TimeEnd));
        set(findobj('tag','etHigh'),'String',num2str(BandpassHigh));
        set(findobj('tag','etLow'),'String',num2str(BandPassLow));
        set(findobj('tag','PSDHz'),'String',num2str(PSDHz));

        set(findobj('tag','customTimeStamp'),'String',num2str(timestamp));
        set(findobj('tag','customTimeStamp'),'Value',timestampOption);

        set(findobj('tag','saveWindowData'),'Value',saveWindowData);
        set(findobj('tag','median'),'Value',medianX);

        set(findobj('tag','waveform'),'Value',WaveForm);
        set(findobj('tag','spectorgram'),'Value',Spectogram);
        set(findobj('tag','writeToFile'),'Value',WriteToFile);
        set(findobj('tag','preserveMemory'),'Value',preserveMemory);
        set(findobj('tag','singleChan'),'Value',singleChan);

        set(findobj('tag','PVL'),'Value',PVL);
        set(findobj('tag','PAL'),'Value',PAL);

        set(findobj('tag','GUIselect'),'Value',GUISelect);

        set(findobj('tag','GainA'),'String',num2str(GainA));
        set(findobj('tag','GainB'),'String',num2str(GainB));
        set(findobj('tag','GainC'),'String',num2str(GainC));
        set(findobj('tag','GainD'),'String',num2str(GainD));

        set(findobj('tag','windowType'),'Value',WindowType);
        set(findobj('tag','windowLength'),'String',num2str(WindowLength));
        set(findobj('tag','overlap'),'String',num2str(Overlap));

        set(findobj('tag','figureSizeX'),'String',num2str(FigureSize(1)));
        set(findobj('tag','figureSizeY'),'String',num2str(FigureSize(2)));
        
        set(findobj('Tag','fontSize'),'String',fontSize);
        set(findobj('Tag','fontType'),'Value',fontType);
        set(findobj('Tag','colormapH'),'Value',colormapH);
        set(findobj('Tag','autoSelectColorbar'),'Value',autoSelectColorbar);
        set(findobj('Tag','publishableFigures'),'Value',publishableFigures);
        
        set(findobj('Tag','colorLp'),'String',colorLp);
        set(findobj('Tag','colorHp'),'String',colorHp);
        set(findobj('Tag','colorLv'),'String',colorLv);
        set(findobj('Tag','colorHv'),'String',colorHv);
        set(findobj('Tag','colorLa'),'String',colorLa);
        set(findobj('Tag','colorHa'),'String',colorHa);
        
        set(findobj('Tag','colorLp'),'Value',AutoSelectFreqRange);
        get(findobj('Tag','freqMin'),'String',freqMin);
        get(findobj('Tag','freqMax'),'String',freqMax);
        
        set(findobj('Tag','consistencyP'),'String',consistencyP);
        set(findobj('Tag','consistencyV'),'String',consistencyV);
        set(findobj('Tag','consistencyA'),'String',consistencyA);
        
        set(findobj('Tag','CalibUnitPressure'),'Value', calibUnitPressure);
        set(findobj('Tag','CalibUnitVelocity'),'Value', calibUnitVelocity);
        
        disp('...Saved Settings loaded...');
%    end

function saveNewProfile(profileName)
    %% Load values from UI
    Threshold = str2double(get(findobj('tag','Threshold'),'String'));
    AutoDetectDelay = str2double(get(findobj('tag','IgnoreFirstSeconds'),'String'));
    FileType = get(findobj('tag','FileType'),'Value');
    AnalysisType = get(findobj('tag','analysisType'),'Value');
    AnalysisParameter = get(findobj('tag','analysisParameter'),'Value');
    SampleRate = str2double(get(findobj('tag','SampleRate'),'String'));
    
    TimeStart = str2double(get(findobj('tag','etTo'),'String'));
    TimeEnd = str2double(get(findobj('tag','etFrom'),'String'));
    BandpassHigh = str2double(get(findobj('tag','etHigh'),'String'));
    BandPassLow = str2double(get(findobj('tag','etLow'),'String'));
    PSDHz = str2double(get(findobj('tag','PSDHz'),'String'));

    saveWindowData = get(findobj('tag','saveWindowData'),'Value');
    medianX = get(findobj('tag','median'),'Value');
    
    timestamp = str2double(get(findobj('tag','customTimeStamp'),'String'));
    timestampOption = get(findobj('tag','customTimeStamp'),'Value');
    
    Stats = 0;
    WaveForm = get(findobj('tag','waveform'),'Value');
    Spectogram = get(findobj('tag','spectorgram'),'Value');
    WriteToFile = get(findobj('tag','writeToFile'),'Value');
    preserveMemory = get(findobj('tag','preserveMemory'),'Value');
    singleChan = get(findobj('tag','singleChan'),'Value');

    PVL = get(findobj('tag','PVL'),'Value');
    PAL = get(findobj('tag','PAL'),'Value');
    percentiles = get(findobj('tag','median'),'Value');
    
    GUISelect = get(findobj('tag','GUIselect'),'Value');
    
    GainA = str2double(get(findobj('tag','GainA'),'String'));
    GainB = str2double(get(findobj('tag','GainB'),'String'));
    GainC = str2double(get(findobj('tag','GainC'),'String'));
    GainD = str2double(get(findobj('tag','GainD'),'String'));
    
    WindowType = get(findobj('tag','windowType'),'Value');
    WindowLength = str2double(get(findobj('tag','windowLength'),'String'));
    Overlap = str2double(get(findobj('tag','overlap'),'String'));
    
    FigureSize(1) = str2double(get(findobj('tag','figureSizeX'),'String'));
    FigureSize(2) = str2double(get(findobj('tag','figureSizeY'),'String'));
    
    fontSize = get(findobj('Tag','fontSize'),'String');
    fontType = get(findobj('Tag','fontType'),'Value');
    colormapH = get(findobj('Tag','colormapH'),'Value');
    autoSelectColorbar = get(findobj('Tag','autoSelectColorbar'),'Value');
    publishableFigures = get(findobj('Tag','publishableFigures'),'Value');

    AutoSelectFreqRange = get(findobj('Tag','colorLp'),'Value');
    freqMin = get(findobj('Tag','freqMin'),'String');
    freqMax = get(findobj('Tag','freqMax'),'String');
    
    colorLp = get(findobj('Tag','colorLp'),'String');
    colorHp = get(findobj('Tag','colorHp'),'String');
    colorLv = get(findobj('Tag','colorLv'),'String');
    colorHv = get(findobj('Tag','colorHv'),'String');
    colorLa = get(findobj('Tag','colorLa'),'String');
    colorHa = get(findobj('Tag','colorHa'),'String');

    consistencyP = get(findobj('Tag','consistencyP'),'String');
    consistencyV = get(findobj('Tag','consistencyV'),'String');
    consistencyA = get(findobj('Tag','consistencyA'),'String');
    
    calibUnitPressure = get(findobj('Tag','CalibUnitPressure'),'Value');
    calibUnitVelocity = get(findobj('Tag','CalibUnitVelocity'),'Value');
    
    try
       load('pathCalib.mat');
    catch
       pathCalib = '';
       calibFile = '';
    end
    
    warning('off','all');
        mkdir('Profiles');
    warning('on','all');
    
    %% Save to file
    if ispc
        profileLocation = ['Profiles\' profileName '.mat'];
    else
        profileLocation = ['Profiles/' profileName '.mat'];
    end
    
    save(profileLocation);
    disp('...settings saved...');


% --- Executes on button press in AutoSelectFreqRange.
function AutoSelectFreqRange_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
       set(findobj('Tag','freqMin'),'Enable','off');
       set(findobj('Tag','freqMax'),'Enable','off');
    else
       set(findobj('Tag','freqMin'),'Enable','on');
       set(findobj('Tag','freqMax'),'Enable','on');
    end

function freqMin_Callback(hObject, eventdata, handles)
% hObject    handle to freqMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqMin as text
%        str2double(get(hObject,'String')) returns contents of freqMin as a double


% --- Executes during object creation, after setting all properties.
function freqMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqMax_Callback(hObject, eventdata, handles)
% hObject    handle to freqMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqMax as text
%        str2double(get(hObject,'String')) returns contents of freqMax as a double


% --- Executes during object creation, after setting all properties.
function freqMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function createMultiChanWavs_Callback(hObject, eventdata, handles)
PM_Analysis_CombineWavFiles

% --------------------------------------------------------------------
function createCalibrationTone_Callback(hObject, eventdata, handles)
PM_Analysis_CreateCalibTone



function ThresholdWait_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdWait (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdWait as text
%        str2double(get(hObject,'String')) returns contents of ThresholdWait as a double


% --- Executes during object creation, after setting all properties.
function ThresholdWait_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdWait (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CalibUnitPressure.
function CalibUnitPressure_Callback(hObject, eventdata, handles)
% hObject    handle to CalibUnitPressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CalibUnitPressure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CalibUnitPressure


% --- Executes during object creation, after setting all properties.
function CalibUnitPressure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalibUnitPressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CalibUnitVelocity.
function CalibUnitVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to CalibUnitVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CalibUnitVelocity contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CalibUnitVelocity


% --- Executes during object creation, after setting all properties.
function CalibUnitVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalibUnitVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
