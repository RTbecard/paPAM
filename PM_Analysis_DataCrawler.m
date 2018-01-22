%Loop current folder and sub folders for .csv files
function PM_Analysis_DataCrawler(filePath,CalibPath,AnalysisType,...
    AnalysisParameter,Threshold,ThresholdWait,GUISelect,...
    Time,BandPass,Output,Gain,SampleRate,timeStamp,windowS,...
    figureOptions,consistencyThresholds,FFTfigure)

    %Consistency analysis parameters
    consistencyP = consistencyThresholds(1);
    consistencyV = consistencyThresholds(2);
    consistencyA = consistencyThresholds(3);

    %Figure settings
    figureS = figureOptions{1};
    fontType = figureOptions{3};
    fontSize = figureOptions{2};
    publishableFigures = figureOptions{4};
    colormapVal = figureOptions{5};
    
    % Colorbar limits
    autoSelectColorbar = figureOptions{6};
    colorLp =  figureOptions{7}(1,1);
    colorHp =  figureOptions{7}(1,2);
    colorLv =  figureOptions{7}(2,1);
    colorHv =  figureOptions{7}(2,2);
    colorLa =  figureOptions{7}(3,1);
    colorHa =  figureOptions{7}(3,2);
    
    % For batch processing
    firstfile = true;
    fileNumber = 0;
    
    % Output options
    outputWavform = Output{1};
    outputSpectrogram = Output{2};
    writeOutputToFile = Output{3};
    PVL = Output{4};
    PAL = Output{5};
    
    if ispc
        outputfolder = [Output{6} '\'];
    else
        outputfolder = [Output{6} '/'];
    end
    singleChannelAnalysis = Output{7};
    outputPercentiles = Output{8};
    
    exportWindowSegments = Output{9};
    preserveMemory = Output{10};
    
    warning('off','all')
    mkdir([outputfolder 'Figures']);
    mkdir([outputfolder 'Results']);
    warning('on','all')
    
    %% Define timeStamp
    if timeStamp{1} == 2;
        timeStamp = timeStamp{2};
    else
       %Date
       tempclock = clock;
       temp = num2str(tempclock(1));
       
       timeStamp = temp;
       temp = num2str(tempclock(2));
       if length(temp)<2
           temp = ['0' temp];
       end
       timeStamp = [timeStamp temp];
       
       temp = num2str(tempclock(3));
       if length(temp)<2
           temp = ['0' temp];
       end
       timeStamp = [timeStamp temp '-'];
       
       %Time
       temp = num2str(tempclock(4));
       if length(temp)<2
           temp = ['0' temp];
       end
       timeStamp = [timeStamp temp];
       
       temp = num2str(tempclock(5));
       if length(temp)<2
           temp = ['0' temp];
       end
       timeStamp = [timeStamp temp];
       
       temp = num2str(int16(tempclock(6)));
       if length(temp)<2
           temp = ['0' temp];
       end
       timeStamp = [timeStamp temp];
    end
    
    %set size for printed figures
    dimentions = figureS;
    
%%  Scan files in current folder
    if ispc()
        if strcmp(filePath{1}(end),'\') == 0;
           filePath{1} =  [filePath{1} '\'];
        end
    else
        if strcmp(filePath{1}(end),'/') == 0;
           filePath{1} =  [filePath{1} '/'];
        end
    end
    
    Filelist = dir([filePath{1} filePath{2}]);
    
    %Filter out folder paths
    tempFilelist = [];
    for l=1:length(Filelist)
        if length(Filelist(l).name) > 4
            if isempty(tempFilelist)
                tempFilelist = Filelist(l);
            else
                tempFilelist(end+1) = Filelist(l);
            end
        end
    end
    
    Filelist = tempFilelist;
    NumberOfFiles = length(Filelist);
    timer = clock;
    
    try
        for i=1:NumberOfFiles
            %% Loop all files in processing list
            FileName = Filelist(i).name;

            %%  Quickly scan file for number of channels
            disp('');
            disp('***************************');
            disp(['File name:' FileName]);
            channels = 0;
            skip3 = false;
            temp = Filelist(i).name;
            if strcmp(temp(length(Filelist(i).name)-3:end),'.csv') ||...
                strcmp(temp(length(Filelist(i).name)-3:end),'.CSV');
               extention = '.csv';
               if ispc()
                    fid = fopen([filePath{1} Filelist(i).name]);
               else
                   fid = fopen([filePath{1} Filelist(i).name]);
               end
               temp = textscan(fid,'%s',1);
               temp = textscan(temp{1}{1},'%s','Delimiter',',');
               channels = size(temp{1},1);
               disp('...csv file selected');
               disp(['Number of Channels:' num2str(channels)]);
            elseif strcmp(temp(length(Filelist(i).name)-3:end),'.wav') || ...
                    strcmp(temp(length(Filelist(i).name)-3:end),'.WAV');
               extention = '.wav';
               % Read the number of channels
               if verLessThan('matlab', '8.1');
                   tempinfo = wavread([filePath{1} Filelist(i).name],1);
               else
                   tempinfo = audioread([filePath{1} Filelist(i).name],[1 1]);
               end
               channels = size(tempinfo,2);
               disp('...wav file selected');
               disp(['Number of Channels:' num2str(channels)]);
            elseif strcmp(temp(length(Filelist(i).name)-3:end),'.mp3') ||...
                    strcmp(temp(length(Filelist(i).name)-3:end),'.MP3');
               extention = '.wav';
               % Read the number of channels
               if verLessThan('matlab', '8.1');
                   disp('Matlab 8.1 or higher required to read .mp3 files!');
                   channels = 0;
                   disp('mp3 file selected...');
                   disp('...skipping file');
               else
                   tempinfo = audioread([filePath{1} Filelist(i).name],[1 1]);
                   channels = size(tempinfo,2);
                   disp('...mp3 file selected');
                   disp(['Number of Channels:' num2str(channels)]);
               end
            else
                disp('Unknown file type :(')
                extention = '';
                channels = 0;
                skip3 = true;
            end

            skip2 = false;
            if (channels < 1) || (channels > 4)
                disp('Error:Incorrect number of channels in file');
                disp('Skipping file');
                skip2 = true;
            end          

            if (PAL == 0) && (PVL == 0) && (channels > 1)
                skip2 = true;
                disp('More than 1 channel, skipping file.');
            end
            % Check for proper filetype
            if skip3 == false && (skip2 == false)
                fileNumber = fileNumber + 1;
                %  Set a timer for file read
                c = clock;
                %% Create filehandles for saving acoustic data to memory
                disp(strcat('Reading file... (File #:',num2str(fileNumber),')'));

                %Velocity
                fidWriteA = fopen('TempCSVDataA.dat','w');
                fidWriteB = fopen('TempCSVDataB.dat','w');
                fidWriteC = fopen('TempCSVDataC.dat','w');

                %Pressure
                fidWriteD = fopen('TempCSVDataD.dat','w');

                if fidWriteA == -1 || fidWriteB == -1 || fidWriteC == -1 || fidWriteD == -1
                    disp('File create error');
                end

                % Open CSV file for reading
                fidRead = fopen([filePath{1} FileName],'r');

                % Convert recorder gain to linear units
                conversionFactorA = Gain(1);
                conversionFactorB = Gain(2);
                conversionFactorC = Gain(3);
                conversionFactorD = Gain(4);

               %Change to hydrophone if specified in GUI
                if channels == 1 && singleChannelAnalysis == 2
                    channels = -1;
                    PVL = 0;
                    PAL = 0;
                end

                %% Copy voltage data into memmap variables (and apply recorder gain corrections)
                %  Using memmap variables stores out data onto the harddrive, so our analysis files are not limited by the about of ram we have.
                %  It takes a bit longer to set up and proccess as a
                %  side-effect
                disp('Copying raw data into temporary space on harddrive...');
                disp('Data will be converted into absolute Voltage.  (Gain factors will be applied to raw data)')
                switch extention
                    case '.csv'
                        % Line-by-line, scan the csv data and copy it to
                        % sererate .dat files on harddrive
                        while ~feof(fidRead)
                            %Convert to voltage units while writing to dat file (Picoscope Oscilloscope) 
                            % These files will be the memmap variables
                            temp = textscan(fidRead, '%f %f %f %f',1e7, 'delimiter',',');
                            % Velocity
                            if channels > 0
                                fwrite(fidWriteA,temp{1,1}.*conversionFactorA,'double');
                            end
                            if channels > 1
                                fwrite(fidWriteB,temp{1,2}.*conversionFactorB,'double');
                            end
                            if channels > 2
                                fwrite(fidWriteC,temp{1,3}.*conversionFactorC,'double');
                            end
                            % Presssure
                            if channels == 4 || channels == -1;
                                fwrite(fidWriteD,temp{1,4}.*conversionFactorD,'double');
                            end
                        end
                    case '.wav'
                        % Get wav file size (needed to know how much hardd drive space we need to store data!)
                        if verLessThan('matlab', '8.1');
                            [temp,wavInfo] = wavfinfo([filePath{1} FileName]); %#ok<*REMFF1>
                            samples = str2double(wavInfo(29:length(wavInfo)-24));
                        else
                            wavInfo = audioinfo([filePath{1} FileName]);
                            samples = int32(wavInfo.TotalSamples); 
                        end
                        samplesIndex = 1;
                        bufferSize = 400000; 
                        while (samplesIndex<samples)
                            %Convert to voltage units while writing to dat file (Picoscope Oscilloscope) 
                            if bufferSize + samplesIndex - 1 > samples
                                linesToWrite = int32(samples - samplesIndex);
                            else
                                linesToWrite = int32(bufferSize);
                            end

                            if verLessThan('matlab', '8.1');
                                [temp,SampleRate] = wavread([filePath{1} FileName], [samplesIndex (samplesIndex + linesToWrite)]);
                            else
                                [temp,SampleRate] = audioread([filePath{1} FileName], double([samplesIndex (samplesIndex + linesToWrite)]));
                            end

                            if channels > 0
                                fwrite(fidWriteA,temp(:,1).*conversionFactorA,'double');
                            end
                            if channels > 1
                                fwrite(fidWriteB,temp(:,2).*conversionFactorB,'double');
                            end
                            if channels > 2
                                fwrite(fidWriteC,temp(:,3).*conversionFactorC,'double');
                            end
                            % Pressure
                            if channels == 4
                                fwrite(fidWriteD,temp(:,4).*conversionFactorD,'double');
                            elseif channels == -1
                                fwrite(fidWriteD,temp(:,1).*conversionFactorD,'double');
                            end
                            samplesIndex = int32(samplesIndex + linesToWrite);
                        end
                end
                %Close the files and delete unecessary variables after
                %the copying/converting is finished.
                clear temp;
                fclose(fidRead);
                fclose(fidWriteA);
                fclose(fidWriteB);
                fclose(fidWriteC);
                fclose(fidWriteD);

                clear fidWriteA fidWriteB fidWriteC fidWriteD fidRead ChanA ChanB ChanC ChanD;

                c2 = clock - c;
                disp(['..file copy time: ',num2str(c2(4)*3600 + c2(5)*60 + fix(c2(6))), ' seconds']);

                % Contaniers for velocity
                if channels > 0;
                    RawDataA = memmapfile('TempCSVDataA.dat','Writable',true,'Format','double');
                end
                if channels > 1
                    RawDataB = memmapfile('TempCSVDataB.dat','Writable',true,'Format','double');
                end
                if channels > 2
                    RawDataC = memmapfile('TempCSVDataC.dat','Writable',true,'Format','double');
                end       
                % Contaniers for pressure
                if channels == -1 || channels == 4
                    RawDataD = memmapfile('TempCSVDataD.dat','Writable',true,'Format','double');
                end

            %%  Start Device calibration
                disp('Creating empty containers for calibrated data...')

            %%  Create empty memmap variables for pressure velocity and acceleration data
                %% Velocity
                if PVL == 1
                    fwriteTemp = fopen('TempCSVDataAbsA.dat','w');
                    fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                    fclose(fwriteTemp);
                    if channels > 1
                        fwriteTemp = fopen('TempCSVDataAbsB.dat','w');
                        fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                        fclose(fwriteTemp);
                    end
                    if channels > 2
                        fwriteTemp = fopen('TempCSVDataAbsC.dat','w');
                        fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                        fclose(fwriteTemp);
                    end
                end
                %% Acceleration
                if channels > 0
                    fwriteTemp = fopen('TempCSVDataAbsAacc.dat','w');
                    fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                    fclose(fwriteTemp);
                end
                if channels > 1
                    fwriteTemp = fopen('TempCSVDataAbsBacc.dat','w');
                    fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                    fclose(fwriteTemp);
                end
                if channels > 2
                    fwriteTemp = fopen('TempCSVDataAbsCacc.dat','w');
                    fwrite(fwriteTemp,zeros(length(RawDataA.Data),1),'double');
                    fclose(fwriteTemp);
                end
                %% Pressure
                if channels > 3 || channels == -1;
                    fwriteTemp = fopen('TempCSVDataAbsD.dat','w');
                    fwrite(fwriteTemp,zeros(length(RawDataD.Data),1),'double');
                    fclose(fwriteTemp);
                end

                clear fwriteTemp;
                % Velocity
                if PVL == 1;
                    AbsDataA = memmapfile('TempCSVDataAbsA.dat','Writable',true,'Format','double');
                    if channels > 1
                        AbsDataB = memmapfile('TempCSVDataAbsB.dat','Writable',true,'Format','double');
                    end
                    if channels > 2
                        AbsDataC = memmapfile('TempCSVDataAbsC.dat','Writable',true,'Format','double');
                    end
                end
                % Acceleration
                if channels > 0
                    AbsDataAacc = memmapfile('TempCSVDataAbsAacc.dat','Writable',true,'Format','double');
                end
                if channels > 1
                    AbsDataBacc = memmapfile('TempCSVDataAbsBacc.dat','Writable',true,'Format','double');
                end
                if channels > 2
                    AbsDataCacc = memmapfile('TempCSVDataAbsCacc.dat','Writable',true,'Format','double');
                end
                % Pressure
                if channels > 3 || channels == -1
                    AbsDataD = memmapfile('TempCSVDataAbsD.dat','Writable',true,'Format','double');
                end

                disp('Copying calibrated data into containers...')

                %% Load calibration file
                % load file
                if ispc;
                    DeviceCalibrationFile = csvread([CalibPath{1} '\' CalibPath{2}]);
                else
                    DeviceCalibrationFile = csvread([CalibPath{1} '/' CalibPath{2}]);
                end
                % remove blank values from file
                for k = 1:4
                    [~,I] = max(DeviceCalibrationFile(:,((k-1)*2) + 1));
                    calibration{k} = DeviceCalibrationFile(1:I,((k-1)*2) + 1:((k-1)*2) + 2);
                end
                % find min frequencies
                switch channels
                    case 1
                        minFrequency = min(min(calibration{1}(:,1)));
                    case 2     
                        minFrequency  = min([min(calibration{1}(:,1)) min(calibration{2}(:,1))]);
                    case 3
                        minFrequency  = min([min(calibration{1}(:,1)) min(calibration{2}(:,1)) min(calibration{3}(:,1))]);
                    case 4
                        minFrequency  = min([min(calibration{1}(:,1)) min(calibration{2}(:,1)) min(calibration{3}(:,1)) min(calibration{4}(:,1))]);
                    case -1
                        minFrequency  = min(min(calibration{4}(:,1)));
                end

                if minFrequency < 1
                    msgbox('Error:  The calibration file has an invalid frequency range.  All frequencies listed in the calibration csv (columns 1,3,5,7) file must be > 0 Hz');
                    disp('Error: Bad calibration file!');
                    return;
                end

                %% Create list of channels to be calibrated
                if channels == -1
                    chanLoop = 4;
                else
                   chanLoop = 1:1:channels;
                end

                %% Open figure for calibration graph
                setCurrentFigure(20);
                hold on;
                
                interpCalibdB = {[0 0] [0 0] [0 0] [0 0]};
                
                %% Loop channels for calibration
                for e=chanLoop
                    %% Define bandpass filter for the calibrated range
                    %first frequency value from each channel
                    switch channels
                        case 1
                            maxCalibrationVal = min([max(calibration{1}(:,1)) (SampleRate./2)-1]);
                            minCalibrationVal = max(calibration{1}(1,1));
                        case 2     
                            maxCalibrationVal = min([max(calibration{1}(:,1)) max(calibration{2}(:,1)) (SampleRate./2)-1]);
                            minCalibrationVal = max([calibration{1}(1,1) calibration{2}(1,1)]);
                        case 3
                            maxCalibrationVal = min([max(calibration{1}(:,1)) max(calibration{2}(:,1)) max(calibration{3}(:,1)) (SampleRate./2)-1]);
                            minCalibrationVal = max([calibration{1}(1,1) calibration{2}(1,1) calibration{3}(1,1)]);
                        case 4
                            maxCalibrationVal = min([max(calibration{1}(:,1)) max(calibration{2}(:,1)) max(calibration{3}(:,1)) max(calibration{4}(:,1)) (SampleRate./2)-1]);
                            minCalibrationVal = max([calibration{1}(1,1) calibration{2}(1,1) calibration{3}(1,1) calibration{4}(1,1)]);
                        case -1
                            maxCalibrationVal = min([max(calibration{4}(:,1)) (SampleRate./2)-1]);
                            minCalibrationVal = max(calibration{4}(1,1));
                    end
                    % Max calibration value is either the max frequency
                    % listed in calibration file or the nyquest frequency

                    % define filter for calibrated range
                    [B,A] = butter(3,[minCalibrationVal maxCalibrationVal]./(SampleRate./2));

                    % Initiate loop
                    samplesIndex = 1;
                    bufferSize = 800000;

                    if channels == -1
                        samples = length(AbsDataD.data);
                    else
                        samples = length(AbsDataAacc.data);
                    end
                 %%  Incrementally calibrate portions of file in loop
                    while (samplesIndex<samples)
                        % Define portion of file to read
                        if bufferSize + samplesIndex - 1 > samples
                            linesToWrite = int32(samples - samplesIndex);
                        else
                            linesToWrite = int32(bufferSize)-1;
                        end
                    %% FFT (turn signal into frequency domain to apply frequency specific calibration constants)
                        switch e;
                            case 1;
                                % Trim calibration file to channel of interest
                                cConstants = calibration{1};
                                % Filter out empty values in calibration file
                                cConstantsIndex = cConstants(:,1) > 0;
                                cConstants = cConstants(cConstantsIndex,:);
                                % run fft on a block of raw data
                                [f0s,Spect]=fftFunc(RawDataA.data(samplesIndex:(samplesIndex + linesToWrite),1),1/SampleRate);
                            case 2;
                                cConstants = calibration{2};
                                cConstantsIndex = cConstants(:,1) > 0;
                                cConstants = cConstants(cConstantsIndex,:);
                                [f0s,Spect]=fftFunc(RawDataB.data(samplesIndex:(samplesIndex + linesToWrite),1),1/SampleRate);
                            case 3;
                                cConstants = calibration{3};
                                cConstantsIndex = cConstants(:,1) > 0;
                                cConstants = cConstants(cConstantsIndex,:);
                                [f0s,Spect]=fftFunc(RawDataC.data(samplesIndex:(samplesIndex + linesToWrite),1),1/SampleRate);
                            case 4;
                                cConstants = calibration{4};
                                cConstantsIndex = cConstants(:,1) > 0;
                                cConstants = cConstants(cConstantsIndex,:);
                                [f0s,Spect]=fftFunc(RawDataD.data(samplesIndex:(samplesIndex + linesToWrite),1),1/SampleRate);
                        end
                        
                        % Equation for converting Omni channel
                        % into SPL
                        % SPL = 20*log10(V)+176
                        % 20*log10(rmsX(RawDataD.Data)) + 176

                    %%  Interpolate calibration values and convert to linear continous units
                        % This will also convet your calibration file to
                        % the save frequency resolution as the fft results
                        % i.e.  You will end up with a unique calibration
                        % constant for each value in the the frequency
                        % spectrum of your raw data.

                        %Convert calibration constants to linear values (from dB ref 1V)
                        linear_cConstants=10.^(cConstants(:,2)./20);
                        interpCalib = spline(cConstants(:,1),linear_cConstants,f0s(1:find(f0s>=maxCalibrationVal,1,'first')));

                        % Write flat line on calibration values outside
                        % calibration range
                        interpCalib((find(f0s>=maxCalibrationVal,1,'first') + 1):length(f0s)) = linear_cConstants(end);
                        interpCalib(1:(find(f0s>=minCalibrationVal,1,'first') + 1)) = linear_cConstants(1);
                        
                        interpCalibdB{e} = [f0s' 20*log10(interpCalib)'];
                        
                        if samplesIndex == 1 && firstfile == true;
                            plotCalibrationValues(20,calibration,interpCalibdB);
                            printFigure(20,['Calibration_' timeStamp]);
                        end
                    %% Convert from V to Pa(or m/s) (apply calibration constants to raw data)
                        %Convert raw data V into uPa (or m/s) for each
                        %frequency
                        if e == 4;
                            % Convert units to uPa
                            Spect = Spect./(interpCalib.');
                        else
                            % Convert units to nm/s and (um/s^2)
                            SpectA = (Spect./(interpCalib.') * 1e6) .* (f0s'.*(2*pi));
                            % Acceleration = w * u
                            % Angular frequency (w) = 2 * pi * frequency
                            Spect = (Spect./(interpCalib.') * 1e9);
                        end
                        % Define length of time domain signal
                        NFFT =(2.*(length(Spect)));
                       % NFFT =(length(Spect));
                    %%  Save converted data and apply bandpass filter to calibration range
                        switch e;
                            case 1;
                                if PVL == 1;
                                    % Convert back to time domain
                                    temp = (1/sqrt(2))*NFFT*ifft(Spect,NFFT,'symmetric');
                                    % Write batch to harddrive
                                    AbsDataA.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                                end
                                % Acceleration
                                temp = (1/sqrt(2))*NFFT*ifft(SpectA,NFFT,'symmetric');
                                AbsDataAacc.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                            case 2;
                                if PVL == 1;
                                    % Convert back to time domain
                                    temp = (1/sqrt(2))*NFFT*ifft(Spect,NFFT,'symmetric');
                                    % Write batch to harddrive
                                    AbsDataB.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                                end
                                % Acceleration
                                temp = (1/sqrt(2))*NFFT*ifft(SpectA,NFFT,'symmetric');
                                AbsDataBacc.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                            case 3;
                                if PVL == 1;
                                    % Convert back to time domain
                                    temp = (1/sqrt(2))*NFFT*ifft(Spect,NFFT,'symmetric');
                                    % Write batch to harddrive
                                    AbsDataC.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                                end
                                % Acceleration
                                temp = (1/sqrt(2))*NFFT*ifft(SpectA,NFFT,'symmetric');
                                AbsDataCacc.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                            case 4;
                                % Convert back to time domain
                                temp = (1/sqrt(2))*NFFT*ifft(Spect,NFFT,'symmetric');
                                % Write batch to harddrive
                                AbsDataD.Data(samplesIndex:(samplesIndex + linesToWrite),1) = filtfilt(B,A,temp(1:linesToWrite+1)); 
                        end
                        samplesIndex = samplesIndex + linesToWrite + 1;
                    end
                end
                
                disp(['Device calibration range:' num2str(minCalibrationVal) '-' num2str(maxCalibrationVal) 'Hz']);
                disp('..Raw data has been converted to absolute units..');
                % Enter analysis details into list
                if firstfile == true;
                    AnalysisList = {{AnalysisParameter,Time,BandPass,GUISelect}};
                end
                sampleNum = 1;

                skip = 0; % this is used for finding the threshold start value
                if (AnalysisType == 2);
            %%  Find beginning of sample (before threshold amplitude is reached)
                    % Batch file analysis
                    % Automatically find start of file using the threshold
                    % Convert threshold into linear nm/(s^2) units
                    linearThreshold = (10^(Threshold/20));
                    disp('..find start of sample from threshold value..');

                % Find length of file
                    if channels == -1
                        fileLength = length(AbsDataD.Data);
                    else
                        fileLength = length(AbsDataAacc.Data);
                    end
                    % Load first 20 seconds of primary track
                    ThreshWait = (ThresholdWait*SampleRate);
                    if ThreshWait >= fileLength
                        error('Threshold wait value is too high.  Exceeds file length.');
                    end
                    idx = (ThreshWait:min(fileLength,ThreshWait + (SampleRate*20)));
                    switch channels
                        case -1
                            temp = AbsDataD.Data(idx,1);
                        case 1
                            temp = AbsDataAacc.Data(idx,1);
                        case 2
                            temp = sqrt(AbsDataAacc.Data(idx,1).^2 + AbsDataBacc.Data(idx,1).^2);
                        case {3,4}
                            temp = sqrt(AbsDataAacc.Data(idx,1).^2 + AbsDataBacc.Data(idx,1).^2 + AbsDataCacc.Data(idx,1).^2);
                    end
                    thresholdDelay = find(abs(temp) >= linearThreshold,1) + ThreshWait;

                    if ~isempty(thresholdDelay)
                        % Threshold found
                        if channels == -1
                            disp(['Threshold uPa detected at: ' num2str(thresholdDelay/SampleRate) ' seconds']);
                        else 
                            disp(['Threshold um/s^2 detected at: ' num2str(thresholdDelay/SampleRate) ' seconds']);
                        end
                    else
                        skip = 1;
                    end
                else
                    % if single file analysis selected, start at first sample
                    thresholdDelay = 0;
                end

                %%  If threshold not found, skip analysis
                if (skip == 1)
                    disp(strcat('Threshold amplitude not reached in first 20 seconds)! Sample:',FileName));
                    if channels == -1
                        disp(['Max uPa in first 10 seconds: ' num2str(num2str(20*log10(max(temp)))) ' dB re 1uPa']);
                    else
                        disp(['Max um/s^2 in first 10 seconds: ' num2str(num2str(20*log10(max(temp)))) ' dB re 1um/s^2']);
                    end
                else
                    loop = true;
                    while loop == true; 
                            %%  Initialize results file
                            %print header (used in all results files)
                            resultsHeader = ['Calibration file:' CalibPath{2} ',' 'Recorder Calibration(X Y Z H):' num2str(conversionFactorA) ' ' num2str(conversionFactorB) ' ' num2str(conversionFactorC) ' ' num2str(conversionFactorD) ',,'];
                            resultsSubHeader = 'File name,Sample Rate,Time range (seconds),Bandwidth (Hz)';
                            %%  Initialize results array
                            results = [FileName ',' num2str(SampleRate)];

                            disp('------------------------------');
                            disp(['Sample ' num2str(sampleNum)]);

                            if (AnalysisList{sampleNum}{4} == 0);
                                %% Single file analysis - Numerically defined start and end times
                                % If times were 0 and 0, use whole track for analysis
                                if (AnalysisList{sampleNum}{2}(2) == 0);
                                    % Add threshold to start time
                                    sampleStart = thresholdDelay + 1;
                                    if channels == -1
                                        sampleEnd = size(AbsDataD.Data,1); 
                                    elseif PVL == 1
                                        sampleEnd = size(AbsDataA.Data,1); 
                                    else
                                        sampleEnd = size(AbsDataAacc.Data,1); 
                                    end
                                else
                                    sampleStart = thresholdDelay + AnalysisList{sampleNum}{2}(1) * SampleRate + 1;
                                    sampleEnd = thresholdDelay + AnalysisList{sampleNum}{2}(2) * SampleRate;
                                end
                            elseif (AnalysisList{sampleNum}{4} == 1);
                                %% Single file analysis - Manually select start and end times (GUI)
                                beep;
                                msgbox('Please select the start of the analysis','Select start','modal');
                                uiwait;
                                setCurrentFigure(18);
                                % show figure from start of threshold to end of file.
                                % User selects start and end using GUI

                                % set x axis
                                if channels == -1
                                    totalDuration = size(RawDataD.Data,1)/SampleRate; 
                                    plotA = PM_Analysis_plotData(AbsDataD.Data(thresholdDelay+1:end),10000);
                                    timeAxis = ((1:size(plotA,1))./size(plotA,1)).*totalDuration;
                                    plot(timeAxis,20*log10(plotA),'b');
                                else
                                    totalDuration = size(RawDataA.Data,1)/SampleRate; 
                                    plotA = PM_Analysis_plotData(AbsDataAacc.Data(thresholdDelay+1:end),10000);
                                    timeAxis = ((1:size(plotA,1))./size(plotA,1)).*totalDuration;
                                    plot(timeAxis,20*log10(plotA),'r');
                                end

                                plot(timeAxis,20*log10(plotA),'r');
                                if channels > 1
                                    plotB = PM_Analysis_plotData(AbsDataBacc.Data(thresholdDelay+1:end),10000);
                                    hold on
                                    plot(timeAxis,20*log10(plotB),'g');
                                    hold off

                                end
                                if channels > 2
                                    plotC = PM_Analysis_plotData(AbsDataCacc.Data(thresholdDelay+1:end),10000);
                                    hold on
                                    plot(timeAxis,20*log10(plotC),'b');
                                    hold off
                                end
                                if channels == -1
                                    ylabel('Amplitude (dB re 1\muPa)');
                                    title('One-sided amplitude for pressure channel');
                                else
                                    ylabel('Amplitude (dB re 1\mum/s^2)');
                                    title('One-sided amplitudes for each channel (X,Y,Z = red, green, blue)');
                                end
                                xlabel('Time (seconds)');

                                setCurrentFigure(18);
                                [x(1),~]  = ginput(1);

                                setCurrentFigure(18);
                                hold on;
                                switch channels
                                    case -1
                                        tempMax = max(max(plotA));
                                        tempMin = min(min(plotA));
                                    case 1
                                        tempMax = max(max(plotA));
                                        tempMin = min(min(plotA));
                                    case 2
                                        tempMax = max(max([plotA plotB]));
                                        tempMin = min(min([plotA plotB]));
                                    case 3
                                        tempMax = max(max([plotA plotB plotC]));
                                        tempMin = min(min([plotA plotB plotC]));
                                    case 4
                                        tempMax = max(max([plotA plotB plotC]));
                                        tempMin = min(min([plotA plotB plotC]));
                                end
                                plot([x(1) x(1)],20.*log10([tempMin tempMax]),':k','LineWidth',2);
                                hold off;
                                beep;
                                msgbox('Please select the end of the analysis','Select end','modal');
                                uiwait;
                                setCurrentFigure(18);
                                [x(2),~]  = ginput(1);
                                sampleStart = floor(x(1)*SampleRate + 1) + thresholdDelay;
                                sampleEnd = floor(x(2)*SampleRate + 1);
                                %Update analysis list (these values will be used on the next file)
                                AnalysisList{sampleNum}{2}(1) = x(1);
                                AnalysisList{sampleNum}{2}(2) = x(2);
                                AnalysisList{sampleNum}{4} = 0;
                                close(18);
                            end     

                            if AnalysisList{sampleNum}{3}(2) == 0
                                AnalysisList{sampleNum}{3}(1) = minCalibrationVal;
                                AnalysisList{sampleNum}{3}(2) = maxCalibrationVal;
                            else
                                AnalysisList{sampleNum}{3}(1) = max([AnalysisList{sampleNum}{3}(1) minCalibrationVal]);
                                AnalysisList{sampleNum}{3}(2) = min([AnalysisList{sampleNum}{3}(2) maxCalibrationVal]);
                            end
                            seconds = ['_' num2str(AnalysisList{sampleNum}{2}(1)) '-' num2str(AnalysisList{sampleNum}{2}(2))];

                            disp(['Time range: ' num2str((sampleStart)/SampleRate) '-' num2str((sampleEnd)/SampleRate) ' seconds']);
                            results = [results ',' num2str(AnalysisList{sampleNum}{2}(1)) '-' num2str(AnalysisList{sampleNum}{2}(2)) ',' num2str(AnalysisList{sampleNum}{3}(1)) '-' num2str(AnalysisList{sampleNum}{3}(2))];
                         %% Extract trimmed selection from file  & apply bandpass filters
                            if sampleEnd > samples;
                               beep;
                               disp('Specified time is outside the range of the recording!');
                               disp(['Track length:' num2str(floor(samples/SampleRate)) ' seconds']);
                               disp('Stopping analysis!!!');
                               return;
                            end

                            % Velocity
                            if PVL == 1;
                                trimDataA = applyBandpass(AbsDataA.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                if channels > 1
                                    trimDataB = applyBandpass(AbsDataB.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                end
                                if channels > 2
                                    trimDataC = applyBandpass(AbsDataC.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                end
                            end
                            % Pressure
                            if channels > 3 || channels == -1;
                                trimDataD = applyBandpass(AbsDataD.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                            end
                            % Acceleration
                            if PAL == 1;
                                if channels > 0
                                    trimDataAacc = applyBandpass(AbsDataAacc.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                end
                                if channels > 1
                                    trimDataBacc = applyBandpass(AbsDataBacc.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                end
                                if channels > 2
                                    trimDataCacc = applyBandpass(AbsDataCacc.Data(sampleStart:sampleEnd),AnalysisList{sampleNum}{3}(1),AnalysisList{sampleNum}{3}(2));
                                end
                            end
                            if AnalysisList{sampleNum}{3}(2) > 0;
                                disp(['Bandpass filter applied(' num2str(AnalysisList{sampleNum}{3}(1)) '-' num2str(AnalysisList{sampleNum}{3}(2)) 'Hz)']);
                            else
                                disp('No Bandpass filter applied');
                            end
                         %% Define window size and overlap for PSD's and Spectrograms
                            if windowS{2} == 0;
                               windowL = SampleRate;
                            else
                               windowL = windowS{2};
                            end
                            switch windowS{1}
                                case 'Hann';
                                    window = hann(windowL);
                                case 'Hamming';
                                    window = hamming(windowL);
                            end
                            noverlap = windowL*(windowS{3}/100);
                            overlap = windowS{3};

                            nfft = windowL;

                            if PAL == 1
                                n = length(trimDataAacc);
                            elseif PVL == 1
                                n = length(trimDataA);
                            elseif channels == -1
                                n = length(trimDataD);
                            end
                            if windowL > n
                                beep;
                               msgbox(['Error:  Your window length (' num2str(windowL) ') is longer than the number of samples in the analysed region (' num2str(n) ').'],'Sample Size','modal');
                               return
                            end
                            if outputWavform == 1;
                                %% Print Waveform (Velocity)
                                %Calculate downsamplerate
                                %Figures will take a long time to load if not
                                %downsampled on large files
                                if channels == -1 || channels == 4;
                                    setCurrentFigure(14);
                                    subplot(1,1,1);
                                    plotWaveform(trimDataD,'p');
                                    title('Pressure H');
                                    printFigure(14,['Waveform(Pressure)_' FileName seconds '_' timeStamp]);
                                end
                                % Find max value for measures
                                if PVL == 1;
                                    switch channels
                                        case 1
                                            maxPlotY = max(max(trimDataA));
                                        case 2
                                            maxPlotY = max(max([trimDataA trimDataB]));
                                        case 3
                                            maxPlotY = max(max([trimDataA trimDataB trimDataC]));
                                        case 4
                                            maxPlotY = max(max([trimDataA trimDataB trimDataC]));
                                    end

                                    setCurrentFigure(2);
                                    switch channels
                                        case 1
                                            subplot(1,1,1);
                                        case 2
                                            subplot(2,1,1);
                                        case 3
                                            subplot(2,2,1);
                                        case 4
                                            subplot(2,2,1);
                                    end
                                    % if loading takes too long, downsample raw data?
                                    % If so, might remove important peaks and give bad visual
                                    % info, also, user will not be able to zoom in.
                                    plotWaveform(trimDataA,'v');
                                    ylim([-maxPlotY maxPlotY]);
                                    title('Velocity X');
                                    if channels > 1
                                        switch channels
                                            case 2
                                                subplot(2,1,2);
                                            case 3
                                                subplot(2,2,2);
                                            case 4
                                                subplot(2,2,2);
                                        end
                                        plotWaveform(trimDataB,'v');
                                        ylim([-maxPlotY maxPlotY]);
                                        title('Velocity Y');
                                    end
                                    if channels > 2
                                        switch channels
                                            case 3
                                                subplot(2,2,3);
                                            case 4
                                                subplot(2,2,3);
                                        end
                                        plotWaveform(trimDataC,'v');
                                        ylim([-maxPlotY maxPlotY]);
                                        title('Velocity Z');  
                                    end 
                                    printFigure(2,['Waveform(Velocity)_' FileName seconds '_' timeStamp]);
                                end

                                % Print acceleration graphs
                                if PAL == 1;
                                    setCurrentFigure(3);
                                     switch channels
                                        case 1
                                            maxPlotY = max(max(trimDataAacc));
                                        case 2
                                            maxPlotY = max(max([trimDataAacc trimDataBacc]));
                                        case 3
                                            maxPlotY = max(max([trimDataAacc trimDataBacc trimDataCacc]));
                                        case 4
                                            maxPlotY = max(max([trimDataAacc trimDataBacc trimDataCacc]));
                                    end

                                    switch channels
                                        case 1
                                            subplot(1,1,1);
                                        case 2
                                            subplot(2,1,1);
                                        case 3
                                            subplot(2,2,1);
                                        case 4
                                            subplot(2,2,1);
                                    end
                                    plotWaveform(trimDataAacc,'a');
                                    ylim([-maxPlotY maxPlotY]);
                                    title('Acceleration X');
                                    if channels > 1
                                        switch channels
                                            case 2
                                                subplot(2,1,2);
                                            case 3
                                                subplot(2,2,2);
                                            case 4
                                                subplot(2,2,2);
                                        end
                                        plotWaveform(trimDataBacc,'a');
                                        ylim([-maxPlotY maxPlotY]);
                                        title('Acceleration Y');
                                    end
                                    if channels > 2
                                        switch channels
                                            case 3
                                                subplot(2,2,3);
                                            case 4
                                                subplot(2,2,3);
                                        end
                                        plotWaveform(trimDataCacc,'a');
                                    ylim([-maxPlotY maxPlotY]);
                                        title('Acceleration Z');  
                                    end
                                    printFigure(3,['Waveform(Acceleration)_' FileName seconds '_' timeStamp]);
                                end
                            end
                            if (outputSpectrogram == 1)
                                %% Spectrogram selected
                                if channels == -1 || channels  == 4
                                    setCurrentFigure(12);
                                    subplot(1,1,1);
                                    drawSpecto({trimDataD},'p');
                                    title('Pressure H');
                                    printFigure(12,['Specgram(Pressure)_' FileName seconds '_Channels_' timeStamp]);
                                end
                                if PVL == 1;
                                    setCurrentFigure(6);
                                    switch channels
                                        case 1
                                            subplot(1,1,1);
                                        case 2
                                            subplot(2,1,1);
                                        case 3
                                            subplot(2,2,1);
                                        case 4
                                            subplot(2,2,1);
                                    end
                                    drawSpecto({trimDataA},'v');
                                    title('Velocity X');
                                    if channels > 1
                                        switch channels
                                            case 2
                                                subplot(2,1,2);
                                            case 3
                                                subplot(2,2,2);
                                            case 4
                                                subplot(2,2,2);
                                        end
                                        drawSpecto({trimDataB},'v');
                                        title('Velocity Y');
                                    end
                                    if channels > 2
                                        switch channels
                                            case 3
                                                subplot(2,2,3);
                                            case 4
                                                subplot(2,2,3);
                                        end
                                        drawSpecto({trimDataC},'v');
                                        title('Velocity Z');
                                    end
                                    % print to file if requested
                                    printFigure(6,['Specgram(Velocity)_' FileName seconds '_Channels_' timeStamp]);
                                    if channels > 1
                                        setCurrentFigure(7);
                                        if channels == 2
                                            subplot(1,1,1);
                                        else
                                            subplot(2,1,1);
                                        end
                                        drawSpecto({trimDataA trimDataB}, 'v');
                                        title('Velocity XY');
                                    end
                                    if channels > 2
                                        subplot(2,1,2);
                                        drawSpecto({trimDataA trimDataB trimDataC}, 'v');
                                        title('Velocity XYZ');
                                    end
                                    if channels > 1
                                        printFigure(7,['Specgram(Velocity)_' FileName seconds '_CombinedChannels_' timeStamp]);
                                    end
                                end
                                if PAL == 1;
                                    % Print acceleration graphs
                                    setCurrentFigure(4);
                                    switch channels
                                        case 1
                                            subplot(1,1,1);
                                        case 2
                                            subplot(2,1,1);
                                        case 3
                                            subplot(2,2,1);
                                        case 4
                                            subplot(2,2,1);
                                    end
                                    drawSpecto({trimDataAacc},'a');
                                    title('Acceleration X');
                                    if channels > 1
                                        switch channels
                                            case 2
                                                subplot(2,1,2);
                                            case 3
                                                subplot(2,2,2);
                                            case 4
                                                subplot(2,2,2);
                                        end
                                        drawSpecto({trimDataBacc},'a');
                                        title('Acceleration Y');
                                    end
                                    if channels > 2
                                        switch channels
                                            case 3
                                                subplot(2,2,3);
                                            case 4
                                                subplot(2,2,3);
                                        end
                                        drawSpecto({trimDataCacc},'a');
                                        title('Acceleration Z');
                                    end
                                    if channels > 1
                                        printFigure(4,['Specgram(Acceleration)1_' FileName seconds '_CombinedChannels_' timeStamp]);
                                    end
                                    % Print combined channels
                                    if channels > 1
                                        setCurrentFigure(5);
                                        if channels == 2
                                            subplot(1,1,1);
                                        else
                                            subplot(2,1,1);
                                        end
                                        drawSpecto({trimDataAacc trimDataBacc},'a');
                                        title('Acceleration XY');
                                        if channels > 2
                                            subplot(2,1,2);
                                            drawSpecto({trimDataAacc trimDataBacc trimDataCacc},'a');
                                            title('Acceleration XYZ');
                                        end
                                        printFigure(5,['Specgram(Acceleration)2_' FileName seconds '_CombinedChannels_' timeStamp]);    
                                    end
                                end
                            end
                            c2 = clock;
                        %%  Analyse subsamples
                            switch AnalysisList{sampleNum}{1};
                                case 1
                                    %% PSD pressure channel
                                    if channels == -1 || channels == 4
                                        disp('PSD Analysis')
                                        setCurrentFigure(13);
                                        subplot(1,1,1);
                                        [pxxD,f,XmedianD,X5D,X95D] = pWelchfunct('D');
                                        psdSettings('p');
                                        title('Pressure H');
                                        printFigure(13,['PSD(Pressure)_' FileName '_' num2str(int32(sampleStart*SampleRate)) '-' num2str(int32(sampleEnd*SampleRate)) '_1_' timeStamp]);
                                        savePSDtoFile({'p'},{'H'},{f},{pxxD},{X5D},{XmedianD},{X95D},'Pressure');
                                    end
                                    %%  PSD analysis (Velocity)
                                    if PVL == 1
                                        if channels > 0;
                                            disp('PSD Analysis')
                                            setCurrentFigure(8);
                                            switch channels
                                                case 1
                                                    subplot(1,1,1);
                                                case 2
                                                    subplot(2,1,1);
                                                case 3
                                                    subplot(2,2,1);
                                                case 4
                                                    subplot(2,2,1);
                                            end
                                            [pxxA,f,XmedianA,X5A,X95A] = pWelchfunct('A');
                                            psdSettings('v');
                                            title('Velocity X');
                                        end
                                        if channels > 1
                                            switch channels
                                                case 2
                                                    subplot(2,1,2);
                                                case 3
                                                    subplot(2,2,2);
                                                case 4
                                                    subplot(2,2,2);
                                            end
                                            [pxxB,f,XmedianB,X5B,X95B] = pWelchfunct('B');
                                            psdSettings('v');
                                            title('Velocity Y');
                                        end
                                        if channels > 2
                                            switch channels
                                                case 3
                                                    subplot(2,2,3);
                                                case 4
                                                    subplot(2,2,3);
                                            end
                                            [pxxC,f,XmedianC,X5C,X95C] = pWelchfunct('C');
                                            psdSettings('v');
                                            title('Velocity Z');
                                        end
                                        % Print combined channels
                                        if channels > 1
                                            setCurrentFigure(9);
                                            if channels == 2
                                                subplot(1,1,1);
                                            else
                                                subplot(2,1,1);
                                            end
                                            [pxxXY,f,XmedianXY,X5XY,X95XY] = pWelchfunct('AB');
                                            psdSettings('v');
                                            title('Velocity XY');
                                        end
                                        if channels > 2
                                            subplot(2,1,2);
                                            [pxxXYZ,f,XmedianXYZ,X5XYZ,X95XYZ] = pWelchfunct('ABC');
                                            psdSettings('v');
                                            title('Velocity XYZ');
                                        end
                                        % Print jpegs
                                        switch channels
                                            case 1
                                                savePSDtoFile({'v'},{'X'},{f},{pxxA},{X5A},{XmedianA},{X95A},'Velocity');
                                            case 2
                                                savePSDtoFile({'v' 'v' 'v'},{'X' 'Y' 'XY'},{f},{pxxA pxxB pxxXY},{X5A X5B X5XY},{XmedianA XmedianB XmedianXY},{X95A X95B X95XY},'Velocity');
                                            case {3,4}
                                                savePSDtoFile({'v' 'v' 'v' 'v' 'v'},{'X' 'Y' 'Z' 'XY' 'XYZ'},{f},{pxxA pxxB pxxC pxxXY pxxXYZ},{X5A X5B X5C X5XY X5XYZ},{XmedianA XmedianB XmedianC XmedianXY XmedianXYZ},{X95A X95B X95C X95XY X95XYZ},'Velocity');
                                        end
                                        if channels > 1
                                            printFigure(9,['PSD(Velocity)_' FileName '_' num2str(int32(sampleStart*SampleRate)) '-' num2str(int32(sampleEnd*SampleRate)) '_2_' timeStamp]);
                                        end
                                        printFigure(8,['PSD(Velocity)_' FileName '_' num2str(int32(sampleStart*SampleRate)) '-' num2str(int32(sampleEnd*SampleRate)) '_1_' timeStamp]);
                                    end
                                    %% Print acceleration graphs (PSD)
                                    if PAL == 1;
                                        setCurrentFigure(10);
                                        switch channels
                                            case 1
                                                subplot(1,1,1);
                                            case 2
                                                subplot(2,1,1);
                                            case 3
                                                subplot(2,2,1);
                                            case 4
                                                subplot(2,2,1);
                                        end
                                        [pxxA,f,XmedianA,X5A,X95A] = pWelchfunct('Aa');
                                        psdSettings('a')
                                        title('Acceleration X');
                                        if channels > 1
                                            switch channels
                                                case 2
                                                    subplot(2,1,2);
                                                case 3
                                                    subplot(2,2,2);
                                                case 4
                                                    subplot(2,2,2);
                                            end
                                            [pxxB,f,XmedianB,X5B,X95B] = pWelchfunct('Ba');
                                            psdSettings('a')
                                            title('Acceleration Y');
                                        end
                                        if channels > 2
                                            switch channels
                                                case 3
                                                    subplot(2,2,3);
                                                case 4
                                                    subplot(2,2,3);
                                            end
                                            [pxxC,f,XmedianC,X5C,X95C] = pWelchfunct('Ca');
                                            psdSettings('a')
                                            title('Acceleration Z');
                                        end
                                        %% Combined channels
                                        if channels > 1
                                            setCurrentFigure(11);
                                            if channels == 2
                                                subplot(1,1,1);
                                            else
                                                subplot(2,1,1);
                                            end
                                            [pxxXY,f,XmedianXY,X5XY,X95XY] = pWelchfunct('ABa');
                                            psdSettings('a')
                                            title('Acceleration XY');
                                        end
                                        if channels > 2
                                            subplot(2,1,2);
                                            [pxxXYZ,f,XmedianXYZ,X5XYZ,X95XYZ] = pWelchfunct('ABCa');
                                            psdSettings('a')
                                            title('Acceleration XYZ');
                                        end
                                        switch channels
                                            case 1
                                                savePSDtoFile({'a'},{'X'},{f},{pxxA},{X5A},{XmedianA},{X95A},'Acceleration');
                                            case 2
                                                savePSDtoFile({'a' 'a' 'a'},{'X' 'Y' 'XY'},{f},{pxxA pxxB pxxXY},{X5A X5B X5XY},{XmedianA XmedianB XmedianXY},{X95A X95B X95XY},'Acceleration');
                                            case {3,4}
                                                savePSDtoFile({'a' 'a' 'a' 'a' 'a'},{'X' 'Y' 'Z' 'XY' 'XYZ'},{f},{pxxA pxxB pxxC pxxXY pxxXYZ},{X5A X5B X5C X5XY X5XYZ},{XmedianA XmedianB XmedianC XmedianXY XmedianXYZ},{X95A X95B X95C X95XY X95XYZ},'Acceleration');
                                        end
                                        % Print jpegs
                                        printFigure(10,['PSD(Acceleration)_' FileName '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_1_' timeStamp]);
                                        if channels > 1
                                            printFigure(11,['PSD(Acceleration)_' FileName '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_2_' timeStamp]);
                                        end
                                    end
                                case 2
                                    %% Impulse
                                    %% 0-peak (max amplitude in selected range)
                                    if (channels > 3) || channels == -1
                                        D = 20*log10(max(abs(trimDataD)));
                                        disp('---Pressure---');
                                        disp(['zero-peak H: ' num2str(D) ' dB ref 1uPa']);

                                        resultsHeader = [resultsHeader ',0-to-peak (dB ref 1uPa)']; %#ok<*AGROW>
                                        resultsSubHeader = [resultsSubHeader ',H'];
                                        results = [results ',' num2str(D)];
                                    end

                                    if PVL == 1
                                        disp('---Velocity---');
                                        if channels > 0
                                            A = 20*log10(max(abs(trimDataA)));
                                            disp(['zero-peak X: ' num2str(A) ' dB ref 1nm/s']);
                                        else
                                            A = 0;
                                        end
                                        if channels > 1
                                            B = 20*log10(max(abs(trimDataB)));
                                            disp(['zero-peak Y: ' num2str(B) ' dB ref 1nm/s']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = 20*log10(max(abs(trimDataC)));
                                            disp(['zero-peak Z: ' num2str(C) ' dB ref 1nm/s']);
                                        else
                                            C = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            XY = 20*log10(max(abs(sqrt(trimDataA.^2 + trimDataB.^2))));
                                            disp(['zero-peak XY: ' num2str(XY) ' dB ref 1nm/s']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            XYZ = 20*log10(max(sqrt(abs(trimDataA.^2 + trimDataB.^2 + trimDataC.^2))));
                                            disp(['zero-peak XYZ: ' num2str(XYZ) ' dB ref 1nm/s']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',0-to-peak (dB ref 1nm/s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        if channels > 0
                                            Aacc = 20*log10(max(abs(trimDataAacc)));
                                            disp(['zero-peak X: ' num2str(Aacc) ' dB ref 1um/s^2']);
                                        else
                                            Aacc = 0;
                                        end
                                        if channels > 1
                                            Bacc = 20*log10(max(abs(trimDataBacc)));
                                            disp(['zero-peak Y: ' num2str(Bacc) ' dB ref 1um/s^2']);
                                        else
                                            Bacc = 0;
                                        end
                                        if channels > 2
                                            Cacc = 20*log10(max(abs(trimDataCacc)));
                                            disp(['zero-peak Z: ' num2str(Cacc) ' dB ref 1um/s^2']);
                                        else
                                            Cacc = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            XYacc = 20*log10(max(abs(sqrt(trimDataAacc.^2 + trimDataBacc.^2))));
                                            disp(['zero-peak XY: ' num2str(XYacc) ' dB ref 1um/s^2']);
                                        else
                                            XYacc = 0;
                                        end
                                        if channels > 2
                                            XYZacc = 20*log10(max(sqrt(abs(trimDataAacc.^2 + trimDataBacc.^2 + trimDataCacc.^2))));
                                            disp(['zero-peak XYZ: ' num2str(XYZacc) ' dB ref 1um/s^2']);
                                        else
                                            XYZacc = 0;
                                        end
                                        resultsHeader = [resultsHeader ',0-to-peak (dB ref 1um/s^2),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(Aacc) ',' num2str(Bacc) ',' num2str(Cacc) ',' num2str(XYacc) ',' num2str(XYZacc)];
                                    end
                                    %% 90% energy envelope
                                    % time for 5-95% energy rise (in power)
                                    if channels == -1 || channels > 3
                                        disp('---Pressure---');
                                        [~,D,rD,SELssD] = calcCulmativeEnergy({trimDataD},4,'p');
                                        resultsHeader = [resultsHeader ',90% energy envelope Pressure (ms),Rise time Pressure (ms),SELss (dB re 1uPa^2*s)'];
                                        resultsSubHeader = [resultsSubHeader ',H,H,H'];
                                        results = [results ',' num2str(D) ',' num2str(rD) ',' num2str(SELssD) ];
                                        printFigure(22,['PulseLength(Pressure)_' FileName '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_1_' timeStamp]);
                                    end
                                    %Velocity
                                    if PVL == 1;
                                        disp('---Velocity---');
                                        if channels > 0
                                            [~,A,rA,SELssA] = calcCulmativeEnergy({trimDataA},1,'v');
                                        else
                                            A = 0;
                                            rA = 0;
                                            SELssA = 0;
                                        end
                                        if channels > 1
                                            [~,B,rB,SELssB] = calcCulmativeEnergy({trimDataB},2,'v');
                                            [~,XY,rXY,SELssXY] = calcCulmativeEnergy({trimDataA trimDataB},5,'v');
                                        else
                                            B = 0;
                                            rB = 0;
                                            SELssB = 0;
                                            XY = 0;
                                            rXY = 0;
                                            SELssXY = 0;
                                        end
                                        if channels > 2
                                            [~,C,rC,SELssC] = calcCulmativeEnergy({trimDataC},3,'v');
                                            [~,XYZ,rXYZ,SELssXYZ] = calcCulmativeEnergy({trimDataA trimDataB trimDataC},6,'v');
                                        else
                                            C = 0;
                                            rC = 0;
                                            SELssC = 0;
                                            XYZ = 0;
                                            rXYZ = 0;
                                            SELssXYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',90% energy envelope velocity (ms),,,,,Rise time velocity (ms),,,,,SELss (dB re (1nm/s)^2*s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ,X,Y,Z,XY,XYZ,X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ) ',' num2str(rA) ',' num2str(rB) ',' num2str(rC) ',' num2str(rXY) ',' num2str(rXYZ) ',' num2str(SELssA) ',' num2str(SELssB) ',' num2str(SELssC) ',' num2str(SELssXY) ',' num2str(SELssXYZ)];
                                        printFigure(21,['PulseLength(Velocity)_' FileName '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_1_' timeStamp]);
                                    end
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        if channels > 0
                                            [~,A,rA,SELssA] = calcCulmativeEnergy({trimDataAacc},1,'a');
                                        else
                                            A = 0;
                                            rA = 0;
                                            SELssA = 0;
                                        end
                                        if channels > 1
                                            [~,B,rB,SELssB] = calcCulmativeEnergy({trimDataBacc},2,'a');
                                            [~,XY,rXY,SELssXY] = calcCulmativeEnergy({trimDataAacc trimDataBacc},5,'a');
                                        else
                                            B = 0;
                                            SELssB = 0;
                                            rB = 0;
                                            XY = 0;
                                            rXY = 0;
                                            SELssXY = 0;
                                        end
                                        if channels > 2
                                            [~,C,rC,SELssC] = calcCulmativeEnergy({trimDataCacc},3,'a');
                                            [~,XYZ,rXYZ,SELssXYZ] = calcCulmativeEnergy({trimDataAacc trimDataBacc trimDataCacc},6,'a');
                                        else
                                            C = 0;
                                            rC = 0;
                                            SELssC = 0;
                                            XYZ = 0;
                                            rXYZ = 0;
                                            SELssXYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',90% energy envelope acceleration (ms),,,,,Rise time acceleration (ms),,,,,SELss (dB re (1um/s)^2*s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ,X,Y,Z,XY,XYZ,X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ) ',' num2str(rA) ',' num2str(rB) ',' num2str(rC) ',' num2str(rXY) ',' num2str(rXYZ) ',' num2str(SELssA) ',' num2str(SELssB) ',' num2str(SELssC) ',' num2str(SELssXY) ',' num2str(SELssXYZ)];
                                        printFigure(19,['PulseLength(Acceleration)_' FileName '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_1_' timeStamp]);
                                    end
                                    %% Crest factor
                                    %20log10(peak/rms)
                                    if channels > 3 || channels == -1
                                        disp('---Pressure---');
                                        D = 20*log10(max(abs(trimDataD))/(rmsX(trimDataD)));
                                        disp(['Crest Factor H: ' num2str(D) ' (dB re 1uPa)']);

                                        resultsHeader = [resultsHeader ',Crest Factor (dB re 1uPa)'];
                                        resultsSubHeader = [resultsSubHeader ',H,X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(D)];
                                    end
                                    % Velocity
                                    if PVL == 1
                                        disp('---Velocity---');
                                        if channels > 0
                                            A = 20*log10(max(abs(trimDataA))/(rmsX(trimDataA)));
                                            disp(['Crest Factor X: ' num2str(A) ' (dB re 1nm/s)']);
                                        else
                                            A = 0;
                                        end
                                        if channels > 1
                                            B = 20*log10(max(abs(trimDataB))/(rmsX(trimDataB)));
                                            disp(['Crest Factor Y: ' num2str(B) ' (dB re 1nm/s)']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = 20*log10(max(abs(trimDataC))/(rmsX(trimDataC)));
                                            disp(['Crest Factor Z: ' num2str(C) ' (dB re 1nm/s)']);
                                        else
                                            C = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            tempData = sqrt(trimDataA.^2 + trimDataB.^2);
                                            tempDataRMS = sqrt(rmsX(trimDataA).^2 + rmsX(trimDataB).^2);
                                            XY = 20*log10(max(abs(tempData))/tempDataRMS);
                                            disp(['Crest Factor XY: ' num2str(XY) ' (dB re 1nm/s)']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            tempData = sqrt(trimDataA.^2 + trimDataB.^2 + trimDataC.^2);
                                            tempDataRMS = sqrt(rmsX(trimDataA).^2 + rmsX(trimDataB).^2 + rmsX(trimDataC).^2);
                                            XYZ = 20*log10(max(abs(tempData))/tempDataRMS);
                                            disp(['Crest Factor XYZ: ' num2str(XYZ) ' (dB re 1nm/s)']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',Crest Factor (dB re 1nm/s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    %Acceleration
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        Aacc = 20*log10(max(abs(trimDataAacc))/(rmsX(trimDataAacc)));
                                        disp(['Crest Factor X (acceleration): ' num2str(Aacc) ' (dB re 1um/s^2)']);
                                        if channels > 1
                                            Bacc = 20*log10(max(abs(trimDataBacc))/(rmsX(trimDataBacc)));
                                            disp(['Crest Factor Y (acceleration): ' num2str(Bacc) ' (dB re 1um/s^2)']);
                                        else
                                            Bacc = 0;
                                        end
                                        if channels > 2
                                            Cacc = 20*log10(max(abs(trimDataCacc))/(rmsX(trimDataCacc)));
                                            disp(['Crest Factor Z (acceleration): ' num2str(Cacc) ' (dB re 1um/s^2)']);
                                        else
                                            Cacc = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            tempData = sqrt(trimDataAacc.^2 + trimDataBacc.^2);
                                            tempDataRMS = sqrt(rmsX(trimDataAacc).^2 + rmsX(trimDataBacc).^2);
                                            XYacc = 20*log10(max(abs(tempData))/tempDataRMS);
                                            disp(['Crest Factor XY (acceleration): ' num2str(XYacc) ' (dB re 1um/s^2)']);
                                        else
                                            XYacc = 0;
                                        end
                                        if channels > 2
                                            tempData = sqrt(trimDataAacc.^2 + trimDataBacc.^2 + trimDataCacc.^2);
                                            tempDataRMS = sqrt(rmsX(trimDataAacc).^2 + rmsX(trimDataBacc).^2 + rmsX(trimDataCacc).^2);
                                            XYZacc = 20*log10(max(abs(tempData))/tempDataRMS);
                                            disp(['Crest Factor XYZ (acceleration): ' num2str(XYZacc) ' (dB re 1um/s^2)']);
                                        else
                                            XYZacc = 0;
                                        end
                                        resultsHeader = [resultsHeader ',Crest Factor (dB re 1um/s^2),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(Aacc) ',' num2str(Bacc) ',' num2str(Cacc) ',' num2str(XYacc) ',' num2str(XYZacc)];
                                    end
                                    printResults('Impulse');
                                case 3 %% Broadband
                                    %% RMS
                                    % Pressure
                                    if channels > 3 || channels == -1
                                        disp('---Pressure---');
                                        D = 20*log10(rmsX((trimDataD)));
                                        disp(['SPL H: ' num2str(D) ' dB re 1uPa']);

                                        resultsHeader = [resultsHeader ',SPL (dB re 1uPa)'];
                                        resultsSubHeader = [resultsSubHeader ',H'];
                                        results = [results ',' num2str(D)];
                                    end

                                    if PVL == 1
                                    % Velocity
                                        disp('---Velocity---');
                                        if channels > 0
                                            A = 20*log10(rmsX((trimDataA)));
                                            disp(['PVL X: ' num2str(A) ' dB re 1nm/s']);
                                        end
                                        if channels > 1
                                            B = 20*log10(rmsX((trimDataB)));
                                            disp(['PVL Y: ' num2str(B) ' dB re 1nm/s']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = 20*log10(rmsX((trimDataC)));
                                            disp(['PVL Z: ' num2str(C) ' dB re 1nm/s']);
                                        else
                                            C = 0;
                                        end
                                        %% Calculate combined channels
                                        if channels > 1
                                            XY = 20*log10(sqrt(rmsX(trimDataA)^2 + rmsX(trimDataB)^2));
                                            disp(['PVL XY: ' num2str(XY) ' dB re 1nm/s']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            XYZ = 20*log10(sqrt(rmsX(trimDataA)^2 + rmsX(trimDataB)^2 + rmsX(trimDataC)^2));
                                            disp(['PVL XYZ: ' num2str(XYZ) ' dB re 1nm/s']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',PVL average (dB re 1nm/s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    %Acceleration
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        Aacc = 20*log10(rmsX((trimDataAacc)));
                                        disp(['PAL X: ' num2str(Aacc) ' dB re (1um/s^2)']);
                                        if channels > 1
                                            Bacc = 20*log10(rmsX((trimDataBacc)));
                                            disp(['PAL Y: ' num2str(Bacc) ' dB re (1um/s^2)']);
                                        else
                                            Bacc = 0;
                                        end
                                        if channels > 2
                                            Cacc = 20*log10(rmsX((trimDataCacc)));
                                            disp(['PAL Z: ' num2str(Cacc) ' dB re (1um/s^2)']);
                                        else
                                            Cacc = 0;
                                        end
                                        if channels > 1
                                            XYacc = 20*log10(sqrt(rmsX(trimDataAacc)^2 + rmsX(trimDataBacc)^2));
                                            disp(['PAL XY: ' num2str(XYacc) ' dB re (1um/s^2)']);
                                        else
                                            XYacc = 0;
                                        end
                                        if channels > 2
                                            XYZacc = 20*log10(sqrt(rmsX(trimDataAacc)^2 + rmsX(trimDataBacc)^2 + rmsX(trimDataCacc)));
                                            disp(['PAL XYZ: ' num2str(XYZacc) ' dB re (1um/s^2)']);
                                        else
                                            XYZacc = 0;
                                        end
                                        resultsHeader = [resultsHeader ',PAL (dB re (1um/s^2),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];     
                                        results = [results ',' num2str(Aacc) ',' num2str(Bacc) ',' num2str(Cacc) ',' num2str(XYacc) ',' num2str(XYZacc)];
                                    end
                                    %% SEL
                                    %Definition of SEL: summed squared sound
                                    %pressure over time
                                    if channels > 3 || channels == -1
                                        disp('---Pressure---');
                                        D = 10*log10(sum(trimDataD.^2)/SampleRate);
                                        disp(['SEL H: ' num2str(D) ' dB re 1uPa^2*s']);

                                        resultsHeader = [resultsHeader ',SEL (dB re (1uPa)^2*s)'];
                                        resultsSubHeader = [resultsSubHeader ',H'];
                                        results = [results ',' num2str(D)];
                                    end
                                    if PVL == 1
                                        if channels > 0
                                            disp('---Velocity---');
                                            A = 10*log10(sum(trimDataA.^2)/SampleRate);
                                            disp(['VEL X: ' num2str(A) ' dB re (1nm/s)^2*s']);
                                        else
                                            A = 0;
                                        end
                                        if channels > 1
                                            B = 10*log10(sum(trimDataB.^2)/SampleRate);
                                            disp(['VEL Y: ' num2str(B) ' dB re (1nm/s)^2*s']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = 10*log10(sum(trimDataC.^2)/SampleRate);
                                            disp(['VEL Z: ' num2str(C) ' dB re (1nm/s)^2*s']);
                                        else
                                            C = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            XY = 10*log10(sqrt(10^(A/10).^2 + 10^(B/10).^2));
                                            disp(['VEL XY: ' num2str(XY) ' dB re (1nm/s)^2*s']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            XYZ = 10*log10(sqrt(10^(A/10).^2 + 10^(B/10).^2 + 10^(C/10).^2));
                                            disp(['VEL XYZ: ' num2str(XYZ) ' dB re (1nm/s)^2*s']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',VEL (dB re (1nm/s)^2*s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        Aacc = 10*log10(sum(trimDataAacc.^2)/SampleRate);
                                        disp(['AEL X: ' num2str(Aacc) ' dB re (1um/s^2)^2*s']);
                                        if channels > 1
                                            Bacc = 10*log10(sum(trimDataAacc.^2)/SampleRate);
                                            disp(['AEL Y: ' num2str(Bacc) ' dB re (1um/s^2)^2*s']);
                                        else
                                            Bacc = 0;
                                        end
                                        if channels > 2
                                            Cacc = 10*log10(sum(trimDataCacc.^2)/SampleRate);
                                            disp(['AEL Z: ' num2str(Cacc) ' dB re (1um/s^2)^2*s']);
                                        else
                                            Cacc = 0;
                                        end
                                        % Combined channels
                                        if channels > 1
                                            XYacc = 10*log10(sqrt(10^(Aacc/10).^2 + 10^(Bacc/10).^2));
                                            disp(['AEL XY: ' num2str(XYacc) ' dB re (1um/s^2)^2*s']);
                                        else
                                            XYacc = 0;
                                        end
                                        if channels > 2
                                            XYZacc = 10*log10(sqrt(10^(Aacc/10).^2 + 10^(Bacc/10).^2 + 10^(Cacc/10).^2));
                                            disp(['AEL XYZ: ' num2str(XYZacc) ' dB re (1um/s^2)^2*s']);
                                        else
                                            XYZacc = 0;
                                        end
                                        resultsHeader = [resultsHeader ',AEL (dB re (1um/s^2)^2*s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(Aacc) ',' num2str(Bacc) ',' num2str(Cacc) ',' num2str(XYacc) ',' num2str(XYZacc)];
                                    end
                                    %% Consistency Analysis (percent time over threshold amplitude)
                                    if channels > 3 || channels == -1% Pressure
                                        disp('---Pressure---');
                                        D = consistencyAnalysis(trimDataD,consistencyP);
                                        disp(['Consistency H: ' num2str(D) ' % above ' num2str(consistencyP) ' dB re 1uPa']);

                                        resultsHeader = [resultsHeader ',Consistency (% above ' num2str(consistencyP) ' dB re 1uPa)'];
                                        resultsSubHeader = [resultsSubHeader ',H'];
                                        results = [results ',' num2str(D)];
                                    end
                                    if PVL == 1 % Velocity
                                        disp('---Velocity---');
                                        if channels > 0
                                            A = consistencyAnalysis(trimDataA,consistencyV);
                                            disp(['Consistency X: ' num2str(A) ' % above ' num2str(consistencyV) ' dB re 1nm/s']);
                                        end
                                        if channels > 1
                                            B = consistencyAnalysis(trimDataB,consistencyV);
                                            disp(['Consistency Y: ' num2str(B) ' % above ' num2str(consistencyV) ' dB re 1nm/s']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = consistencyAnalysis(trimDataC,consistencyV);
                                            disp(['Consistency Z: ' num2str(C) ' % above ' num2str(consistencyV) ' dB re 1nm/s']);
                                        else
                                            C = 0;
                                        end
                                        %% Calculate combined channels
                                        if channels > 1
                                            temp = sqrt(trimDataA.^2 + trimDataB.^2);
                                            XY = consistencyAnalysis(temp,consistencyV);
                                            clear temp;
                                            disp(['Consistency XY: ' num2str(XY) ' % above ' num2str(consistencyV) ' dB re 1nm/s)']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            temp = sqrt(trimDataA.^2 + trimDataB.^2 + trimDataC.^2);
                                            XYZ = consistencyAnalysis(temp,consistencyV);
                                            clear temp;
                                            disp(['Consistency XYZ: ' num2str(XYZ) ' % above ' num2str(consistencyV) ' dB re 1nm/s']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',Consistency (% above ' num2str(consistencyV) ' dB re 1nm/s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    %Acceleration
                                    if PAL == 1;
                                        disp('---Acceleration---');
                                        if channels > 0
                                            A = consistencyAnalysis(trimDataAacc,consistencyA);
                                            disp(['Consistency X: ' num2str(A) ' % above ' num2str(consistencyA) ' dB re 1um/s']);
                                        end
                                        if channels > 1
                                            B = consistencyAnalysis(trimDataBacc,consistencyA);
                                            disp(['Consistency Y: ' num2str(B) ' % above ' num2str(consistencyA) ' dB re 1um/s']);
                                        else
                                            B = 0;
                                        end
                                        if channels > 2
                                            C = consistencyAnalysis(trimDataCacc,consistencyA);
                                            disp(['Consistency Z: ' num2str(C) ' % above ' num2str(consistencyA) ' dB re 1um/s']);
                                        else
                                            C = 0;
                                        end
                                        %% Calculate combined channels
                                        if channels > 1
                                            temp = sqrt(trimDataAacc.^2 + trimDataBacc.^2);
                                            XY = sum(temp < -1*(10^(consistencyV/20))); %num of elements less than thresh
                                            XY = XY + sum(temp > (10^(consistencyV/20))); %num of elements greater than thresh
                                            XY = (XY/length(temp))*100; % convert to percent
                                            clear temp;
                                            disp(['Consistency XY: ' num2str(XY) ' % above ' num2str(consistencyV) ' dB re 1um/s']);
                                        else
                                            XY = 0;
                                        end
                                        if channels > 2
                                            temp = sqrt(trimDataAacc.^2 + trimDataBacc.^2 + trimDataCacc.^2);
                                            XYZ = sum(temp < -1*(10^(consistencyV/20))); %num of elements less than thresh
                                            XYZ = XYZ + sum(temp > (10^(consistencyV/20))); %num of elements greater than thresh
                                            XYZ = (XYZ/length(temp))*100; % convert to percent
                                            clear temp;
                                            disp(['Consistency XYZ: ' num2str(XY) ' % above ' num2str(consistencyV) ' dB re 1um/s']);
                                        else
                                            XYZ = 0;
                                        end
                                        resultsHeader = [resultsHeader ',Consistency (% above ' num2str(consistencyP) ' dB re 1um/s),,,,'];
                                        resultsSubHeader = [resultsSubHeader ',X,Y,Z,XY,XYZ'];
                                        results = [results ',' num2str(A) ',' num2str(B) ',' num2str(C) ',' num2str(XY) ',' num2str(XYZ)];
                                    end
                                    printResults('Broadband');
                            end
                            %if at end of analysis list, exit loop
                            if sampleNum == length(AnalysisList);
                                loop = false;
                            else 
                                loop = true;
                            end
                            % The user can only add extra paramaters in the first file
                            % of batch analysis.  All proceeding files will have the
                            % same analysis list used.
                            if firstfile == true && sampleNum == length(AnalysisList);
                                answer = questdlg('Add another analysis to this file?','Add Analysis','Yes','No','No');
                                if strcmp(answer,'Yes');
                                    % call GUI for input for new analysis paramaters
                                    AnalysisParam = addAnalysis;
                                    for s = 1:size(AnalysisParam,1)
                                        AnalysisList{end + 1} = AnalysisParam(s,:);
                                    end
                                    if sampleNum == length(AnalysisList);
                                        loop = false;
                                    else 
                                        loop = true;
                                    end
                                end
                            end
                            sampleNum = sampleNum + 1;
                    end
                %%  Write final timeStamps
                    disp('***************************')
                    disp(['Finished processing file: ', FileName]);
                    c2 = clock - c2;
                    disp(strcat('Analysis time = ',num2str(c2(4)*3600 + c2(5)*60 + fix(c2(6))),' seconds'));
                    c2 = clock - c;
                    disp(strcat('Total file handling time = ',num2str(c2(4)*3600 + c2(5)*60 + fix(c2(6))),' seconds'));
                    clear A B C D trimDataA trimDataB trimDataC trimDataD;
                    c3 = clock - timer;
                    disp(strcat('Total batch analysis time = ',num2str(c3(4)*3600 + c3(5)*60 + fix(c3(6))),' seconds'));

                    disp(['Percent done: ' num2str(floor((i/length(Filelist))*100))]);
                    disp(['Estimated time remaining: ' num2str((length(Filelist)-i) * (c3(4)*3600 + c3(5)*60 + fix(c3(6)))) ' seconds']);
                    disp('------------------------------');
                    if ispc()
                        memory;
                    end
                    disp('------------------------------');
                    firstfile = false;
                    
                    %% Append PSD window data
                    if exportWindowSegments == 1;
                       combinePSDs(outputfolder,FileName,timeStamp); 
                    end
                end
            end
            disp('Deleting temporary files...');
            clear RawDataA RawDataB RawDataC RawDataD...
                AbsDataA AbsDataB AbsDataC AbsDataD...
                AbsDataAacc AbsDataBacc AbsDataCacc AbsDataDacc;
            warning('off','all');
            delete 'TempCSVDataA.dat' 'TempCSVDataB.dat' 'TempCSVDataC.dat' 'TempCSVDataD.dat'...
                'TempCSVDataAacc.dat' 'TempCSVDataBacc.dat' 'TempCSVDataCacc.dat'...
                'TempCSVDataAbsA.dat' 'TempCSVDataAbsB.dat' 'TempCSVDataAbsC.dat' 'TempCSVDataAbsD.dat'...
                'TempCSVDataAbsAacc.dat' 'TempCSVDataAbsBacc.dat' 'TempCSVDataAbsCacc.dat';
            warning('on','all');
        end
    catch ME
        switch ME.identifier
            case 'MATLAB:nomem'
                msgbox('Out of memory error.  Please check the Memory Issues section of the manual for advice on how to work around this error.');
            otherwise
                ME.rethrow;
        end
        disp('!!!Exiting analysis due to error!!!');
    end
    % Merge results in tempPSD folder into a single file
    appendPSDResults;
    % Clear variables and exit analysis
    clear global;
    clear all;
    disp('------------------------------');
    beep;
    
%%  Data plotting nested function  
    function output = PM_Analysis_plotData(data,bins)
        %Breaks the data into a specified number of bins then returns the
        %abs(maximum) value of each bin.  This is for dramatically speeding up
        %the time it takes to plot amplitude functions of given recordings
        
        %Subsample matrix for image
        subsample = fix(size(data,1)./bins);

        output = zeros(bins,1);
        
        for t=1:bins
            output(t) = max(abs(data((t-1)*subsample + 1:(t)*subsample)));
        end
    end
%%  FFT nested function
    function [freq, P]=fftFunc(p,dt)
        % calculates FFT of given signal, used in the calibration
        % p= signal wave form
        % dt= time step

        nsamp = length(p); %length of signal
        NFFT=nsamp;
        %NFFT = 2^nextpow2(nsamp); % #FFT points(Next power of 2 from length)

        P=fft(p,NFFT);  % FFT of signal
        
        %Disable integer warning
        warning('off','all');
        P=sqrt(2)*P(1:NFFT/2+1)/NFFT; % positive frequencies and normalization
     %  P=dt*NFFT*P(1:NFFT/2+1)/nsamp; % positive frequencies and normalization
     %  Reenable all warnings
        warning('on','all');
        freq = (0.5/dt)*linspace(0,1,NFFT/2+1); % frequency vector

        % NFFT = 2^nextpow2(nsamp); % Next power of 2 from length of y
        % Y= fft(p,NFFT)/nsamp;
        % freq = (1/dt)/2*linspace(0,1,NFFT/2);
        % 
        % P=2*abs(Y(1:NFFT/2)); 
    end

    function printFigure(fig,name)
        if writeOutputToFile == 1
            setCurrentFigure(fig);
            if ispc
                path = [outputfolder 'Figures\' name];
            else
                path = [outputfolder 'Figures/' name];
            end 
            figure(fig);
            if publishableFigures == 1;
                myaaFig = myaa;
            end
            saveas(fig,[path '.jpeg']);
            saveas(fig,[path '.fig']);
            if publishableFigures == 1;
                close(myaaFig);
            end
        end
        
        if preserveMemory == 1
            close(fig);
            try
                uiwait(fig);
            end
        end
    end

    function setFigureSize
        temph = gcf;
        tempfig = get(temph,'Position');
        %Resize only if figure is not already correct size
        if tempfig(3) ~= dimentions(1) || tempfig(4) ~= dimentions(2) 
               set(temph,'Position',[0 0 dimentions(1) dimentions(2)]);
        end
    end

    function drawSpecto(cellData,type)
        if length(cellData) == 1
            [power,F,T] = SingleSpec(cellData{1});
            dB = power;
        elseif length(cellData) == 2
            [power1,F,T] = SingleSpec(cellData{1});
            [power2,~,~] = SingleSpec(cellData{2});
            dB = sqrt(power1.^2 + power2.^2);
        elseif length(cellData) == 3
            [power1,F,T] = SingleSpec(cellData{1});
            [power2,~,~] = SingleSpec(cellData{2});
            [power3,~,~] = SingleSpec(cellData{3});
            dB = sqrt(power1.^2 + power2.^2 + power3.^2);
        end
        
        %print spectogram
        spectroH = imagesc(T,F,real(10*log10(dB)));
        set(get(spectroH,'parent'),'ylim',[minCalibrationVal maxCalibrationVal]);
        view(2);
        specSettings(spectroH);
        
        function [power,F,T] = SingleSpec(data)
            %calculate how many parts to divide spectogram into
            segments = fix(length(data)/400000);
            if segments < 1
                segments = 1;
            end
            segmentLength = fix(length(data)/segments);
            %segments of around 400,000 samples

            %Create segments
            spec = {};
            for p=1:segments
                [temp1,temp2,temp3] = specgram(data((segmentLength*(p-1)) + 1:segmentLength*p),nfft,SampleRate,window,noverlap);
                spec{p} = {temp1 temp2 temp3};
            end
            
            %Stitch segments together
            power = [];
            F = [];
            T = [];
            for p=1:segments
                temp = spec{p};
                power = [power temp{1}];
                F = temp{2};
                if isempty(T)
                    T = temp{3};
                else
                    T = [T; (temp{3} + T(end))];
                end
            end
        end
        function specSettings(spectroH)
            axis tight;
            parentH = get(spectroH,'parent');
            xlabel(get(spectroH,'parent'),'Time (s)');
            ylabel(get(spectroH,'parent'),'Frequency (Hz)');
            set(get(spectroH,'parent'),'YDir','normal');
            colormap(get(spectroH,'parent'),colormapVal);
            t = colorbar('peer',get(spectroH,'parent'));
            switch type
                case 'p'
                    set(get(t,'YLabel'),'String','PSD (dB re (1\muPa)^2/Hz)');
                    if autoSelectColorbar == 0
                        set(get(spectroH,'parent'), 'CLim', [colorLp, colorHp]);
                    end
                case 'v'
                    set(get(t,'YLabel'),'String','PSD (dB re (1nm/s)^2/Hz)');
                    if autoSelectColorbar == 0
                        set(get(spectroH,'parent'), 'CLim', [colorLv, colorHv]);
                    end
                case 'a'
                    set(get(t,'YLabel'),'String','PSD (dB re (1\mum/s^2)^2/Hz)');
                    if autoSelectColorbar == 0
                        set(get(spectroH,'parent'), 'CLim', [colorLa, colorHa]);
                    end
            end
            ylim(get(spectroH,'parent'),[AnalysisList{sampleNum}{3}(1) AnalysisList{sampleNum}{3}(2)]);
        end
    end

    function psdSettings(type)
        xlabel('Frequency (Hz)');
        switch type
            case 'p'        
                ylabel('PSD sound pressure level (dB re (1\muPa)^2/Hz)');
            case 'v'
                ylabel('PSD velocity level (dB re (1nm/s)^2/Hz)');
            case 'a'
                ylabel('PSD acceleration level (dB re (1\mum/s^2)^2/Hz)');
        end
        grid on;
    end

    function [output] = addAnalysis
        % Asks user to input new analysis paramaters
        global xAddAnalysis;
        
        % Clear old values
        xAddAnalysis = [];
        
        open 'PM_Analysis_AddParameter.fig'
        uiwait;

        output = xAddAnalysis;
    end

    function output = applyBandpass(datax,min,max)
        if max > 0;
            % if upper frequency specified apply analysis
            try
                [b,a] = butter(3,[min max]./(SampleRate./2));
            catch
               msgbox('You selected bandpass frequencies must be greater than 0 and less then the nyquest frequency (sample rate / 2).  This sample will be skipped');
                return;
            end
            output = filter(b,a,datax);
        else
            % if no bandpass specified, return origional data
            output = datax;
        end
    end

    function savePSDtoFile(type,channel,f,Xmean,X5,X50,X95,Name)
        if writeOutputToFile == 1
            if ispc
                temppath = [outputfolder 'Results\TempPSD\'];
                warning('off','all');
                mkdir(temppath(1:end-1));
                warning('on','all');
            else
                temppath = [outputfolder 'Results/TempPSD/'];
                warning('off','all');
                mkdir(temppath(1:end-1));
                warning('on','all');
            end

            tempCell = f{1};
            tempColumnName = 'Frequency (Hz)';
            tempName = [temppath 'PSD_' FileName '_' Name '_' num2str(int32(AnalysisList{sampleNum}{3}(1))) '-' num2str(int32(AnalysisList{sampleNum}{3}(2))) '_' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate)) '_' timeStamp '.csv'];
            %% Loop channels
            for w = 1:length(type);
                switch type{w}
                    case 'p'
                        tempColumnName = [tempColumnName ',PSD ' channel{w} ' (dB ref 1uPa^2/Hz)'];
                    case 'v'
                        tempColumnName = [tempColumnName ',PSD ' channel{w} ' (dB ref (1nm/s)^2/Hz)'];
                    case 'a'
                        tempColumnName = [tempColumnName ',PSD ' channel{w} ' (dB ref (1um/s^2)^2/Hz)'];
                end
                if outputPercentiles == 1
                    tempColumnName = [tempColumnName ',5th percentile,50th percentile,95th percentile'];
                    tempCell = [tempCell 10*log10(Xmean{w}) 10*log10(X5{w}) 10*log10(X50{w}) 10*log10(X95{w})];
                else
                    tempCell = [tempCell 10*log10(Xmean{w})];
                end
            end

            writerh = fopen(tempName,'w');
            fprintf(writerh,'%s\n',FileName);
            fprintf(writerh,'%s\n',['Time Range (s):' num2str(int32(sampleStart/SampleRate)) '-' num2str(int32(sampleEnd/SampleRate))]);
            fprintf(writerh,'%s\n',['Bandpass (Hz):' num2str(int32(AnalysisList{sampleNum}{3}(1))) '-' num2str(int32(AnalysisList{sampleNum}{3}(2)))]);
            fprintf(writerh,'%s\n',['Calibration file:' CalibPath{2}]);
            fprintf(writerh,'%s\n',['Recorder Calibration(X Y Z H):' num2str(conversionFactorA) ' ' num2str(conversionFactorB) ' ' num2str(conversionFactorC) ' ' num2str(conversionFactorD)]);
            fprintf(writerh,'%s\n',tempColumnName);
            fclose(writerh);

            Imin = find(tempCell(:,1) >= AnalysisList{sampleNum}{3}(1),1);
            Imax = find(tempCell(:,1) >= AnalysisList{sampleNum}{3}(2),1);

            % Print data to file (in dB units)
            dlmwrite(tempName,[tempCell(Imin:Imax,1) tempCell(Imin:Imax,2:end)],'-append');
        end
    end

    function [pxx,freq,Xmedian,X5,X95] = pWelchfunct(channel)
        switch channel
            % if multiple channels, combine powersepctrums from each
            case 'D'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('D',window,overlap,windowL,SampleRate);
            case 'A'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('A',window,overlap,windowL,SampleRate);
            case 'B'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('B',window,overlap,windowL,SampleRate);
            case 'C'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('C',window,overlap,windowL,SampleRate);
            case 'AB'
                [freq,pxx1,Xmedian1,X51,X951] = PM_Analysis_fft('A',window,overlap,windowL,SampleRate);
                [~,pxx2,Xmedian2,X52,X952] = PM_Analysis_fft('B',window,overlap,windowL,SampleRate);
                pxx = sqrt(pxx1.^2 + pxx2.^2);
                Xmedian = sqrt(Xmedian1.^2 + Xmedian2.^2);
                X5 = sqrt(X51.^2 + X52.^2);
                X95 = sqrt(X951.^2 + X952.^2);
            case 'ABC'
                [freq,pxx1,Xmedian1,X51,X951] = PM_Analysis_fft('A',window,overlap,windowL,SampleRate);
                [~,pxx2,Xmedian2,X52,X952] = PM_Analysis_fft('B',window,overlap,windowL,SampleRate);
                [~,pxx3,Xmedian3,X53,X953] = PM_Analysis_fft('C',window,overlap,windowL,SampleRate);
                pxx = sqrt(pxx1.^2 + pxx2.^2 + pxx3.^2);
                Xmedian = sqrt(Xmedian1.^2 + Xmedian2.^2 + Xmedian3.^2);
                X5 = sqrt(X51.^2 + X52.^2 + X53.^2);
                X95 = sqrt(X951.^2 + X952.^2 + X953.^2);
            case 'Aa'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('Aa',window,overlap,windowL,SampleRate);
            case 'Ba'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('Ba',window,overlap,windowL,SampleRate);
            case 'Ca'
                [freq,pxx,Xmedian,X5,X95] = PM_Analysis_fft('Ca',window,overlap,windowL,SampleRate);
            case 'ABa'
                [freq,pxx1,Xmedian1,X51,X951] = PM_Analysis_fft('Aa',window,overlap,windowL,SampleRate);
                [~,pxx2,Xmedian2,X52,X952] = PM_Analysis_fft('Ba',window,overlap,windowL,SampleRate);
                pxx = sqrt(pxx1.^2 + pxx2.^2);
                Xmedian = sqrt(Xmedian1.^2 + Xmedian2.^2);
                X5 = sqrt(X51.^2 + X52.^2);
                X95 = sqrt(X951.^2 + X952.^2);
            case 'ABCa'
                [freq,pxx1,Xmedian1,X51,X951] = PM_Analysis_fft('Aa',window,overlap,windowL,SampleRate);
                [~,pxx2,Xmedian2,X52,X952] = PM_Analysis_fft('Ba',window,overlap,windowL,SampleRate);
                [~,pxx3,Xmedian3,X53,X953] = PM_Analysis_fft('Ca',window,overlap,windowL,SampleRate);
                pxx = sqrt(pxx1.^2 + pxx2.^2 + pxx3.^2);
                Xmedian = sqrt(Xmedian1.^2 + Xmedian2.^2 + Xmedian3.^2);
                X5 = sqrt(X51.^2 + X52.^2 + X53.^2);
                X95 = sqrt(X951.^2 + X952.^2 + X953.^2);
        end
        %% print PSD with legend
        if outputPercentiles == 1
            h2 = semilogx(freq,10*log10(X5),'--r');
            hold on
            semilogx(freq,10*log10(X95),'--r');
            h3 = semilogx(freq,10*log10(Xmedian),'--g');
            h1 = semilogx(freq,10*log10(pxx));
            hold off
            legend([h1,h2,h3],'Mean','5^t^h and 95^t^h percentile','50^t^h percentile','Location','southeast');
        else
            semilogx(freq,10*log10(pxx));
        end
        
        %% Apply fft figure options
        if FFTfigure{1} == 0 && (FFTfigure{2} < FFTfigure{3})
            Fmin = max(FFTfigure{2},minCalibrationVal);
            Fmax = min(FFTfigure{3},maxCalibrationVal);
            xlim([Fmin Fmax]);
            
            iStart = find(freq >= Fmin,1,'first');
            iEnd = find(freq >=Fmax,1,'first');
            ylim([10*log10(min(X5(iStart:iEnd))) 10*log10(max(X95(iStart:iEnd)))]);
        end
    end

    function plotWaveform(data,type)
        totalDuration = length(data)/SampleRate;
        timeAxis = ((1:length(data))./length(data)).*totalDuration;
        
        plot(timeAxis,data(1:end));
        grid on;
        xlabel('Time (Seconds)');
        switch type 
            case 'p'
                ylabel('Amplitude (\muPa)');
            case 'a'
                ylabel('Amplitude (\mum/s^2)');
            case 'v'
                ylabel('Amplitude (nm/s)');
        end
    end
    function [tempEnergy,pulseLength,riseTime,SELss] = calcCulmativeEnergy(dataset,currentChannel,type)
        tempEnergy = {};
        for w = 1:length(dataset)
            data = dataset{w};
            for p = 1:length(data);
                if p == 1
                    tempEnergy{w}(p) = data(p)^2;
                else
                    tempEnergy{w}(p) = tempEnergy{w}(p-1) + data(p)^2;
                end
            end
            tempEnergy{w} = tempEnergy{w}./SampleRate;
        end
        
        switch length(dataset)
            case 1
                Amp = abs(dataset{1});
                tempEnergy = tempEnergy{1};
            case 2
                Amp = sqrt(dataset{1}.^2 + dataset{2}.^2);
                tempEnergy = sqrt(tempEnergy{1}.^2 + tempEnergy{2}.^2);
            case 3
                Amp = sqrt(dataset{1}.^2 + dataset{2}.^2 + dataset{3}.^2);
                tempEnergy = sqrt(tempEnergy{1}.^2 + tempEnergy{2}.^2 + tempEnergy{3}.^2);
        end
        
        maxAmp = max(abs(Amp));
        
        maxEnergy = tempEnergy(end); 
        maxAmplitude = max(data);
        EI(1) = find(tempEnergy >= maxEnergy*0.05,1);
        EI(2) = find(tempEnergy >= maxEnergy*0.95,1);
        AI(1) = find(Amp >= maxAmp*0.05,1);
        AI(2) = find(Amp >= maxAmp*0.95,1);
        
        pulseLength = num2str((EI(2) - EI(1))*1000/SampleRate);
        riseTime = num2str((AI(2) - AI(1))*1000/SampleRate);
        
        %% Calc SEL for each channel
        SELss = {};
        for w = 1:length(dataset)
            data = dataset{w};
            SELss{w} = sum(data(EI(1):EI(2)).^2)/SampleRate;
        end
        %% Join SEL channels
        switch length(dataset)
            case 1
                SELss = SELss{1};
            case 2
                SELss = sqrt(SELss{1}.^2 + SELss{2}.^2);
            case 3
                SELss = sqrt(SELss{1}.^2 + SELss{2}.^2 + SELss{3}.^2);
        end
        % Convert to dB units
        SELss = 10*log10(SELss);
        
        chans = {'X' 'Y' 'Z' 'H' 'XY' 'XYZ'};
        tempstring = [chans{currentChannel} ': ' num2str(SELss)];
        
        %% Calc SEL single strike
        disp(['Pulse length: ' chans{currentChannel} ': ' pulseLength ' ms']);
        disp(['Rise Time: ' chans{currentChannel} ': ' riseTime ' ms']);
        switch type
            case 'p'
                tempstring = ['SELss ' tempstring ' dB re 1uPa^2*s'];
            case 'v'
                tempstring = ['VELss ' tempstring ' dB re (1nm/s)^2*s'];
            case 'a'
                tempstring = ['AELss ' tempstring ' dB re (1um/s^2)^2*s'];
        end
        disp(tempstring); 
        %% load figure
        
        if strcmp(type,'a') == 1
            setCurrentFigure(19);
        elseif strcmp(type,'v') == 1
            setCurrentFigure(21);
        else
            setCurrentFigure(22);
        end
        
        if currentChannel < 5
            %% Single channels
            if currentChannel == 4
                subplot(1,1,1);
                plotCulmativeEnergy;
            elseif channels == 1 || channels == -1
                subplot(1,1,1);
                plotCulmativeEnergy;
            elseif channels == 2
                subplot(2,1,currentChannel);
                plotCulmativeEnergy;
            elseif channels == 3
                subplot(2,2,currentChannel);
                plotCulmativeEnergy;
            elseif channels == 4
                subplot(2,2,currentChannel);
                plotCulmativeEnergy;
            end
        end

        function plotCulmativeEnergy
            [hAx,~,H2] = plotyy((1:length(data))/SampleRate,data,(1:length(data))/SampleRate,tempEnergy);
            xlim(hAx(2),[0 length(data)/SampleRate]);
            xlim(hAx(1),[0 length(data)/SampleRate]);
            ylim(hAx(1),[-maxAmplitude maxAmplitude]);
            ylim(hAx(2),[0 maxEnergy]);
            hold on
            plot([EI(1)/SampleRate EI(1)/SampleRate],[-maxAmplitude,maxAmplitude],'r');
            plot([AI(1)/SampleRate AI(1)/SampleRate],[-maxAmplitude,maxAmplitude],'m');
            plot([EI(2)/SampleRate EI(2)/SampleRate],[-maxAmplitude,maxAmplitude],'r');
            plot([AI(2)/SampleRate AI(2)/SampleRate],[-maxAmplitude,maxAmplitude],'m');
            xlim([0 length(data)/SampleRate]);
            setLabels;
            hold off
            function setLabels
                title([chans{currentChannel} ' Red = Pulse length; Magenta = Rise Time']);
                xlabel('Seconds');
                set(H2,'linewidth',2);
                switch type 
                    case 'p'
                        ylabel(hAx(1),'Amplitude (\muPa)');
                        ylabel(hAx(2),'Cumulative Power (\muPa^2)*s');
                    case 'a'
                        ylabel(hAx(1),'Amplitude (\mum/s^2)');
                        ylabel(hAx(2),'Cumulative Power (\mum/s^2)^2*s)');
                    case 'v'
                        ylabel(hAx(1),'Amplitude (nm/s)');
                        ylabel(hAx(2),'Cumulative Power (nm/s)^2*s');
                end
            end
        end
    end

    function appendPSDResults
        if ispc
            temppath = [outputfolder 'Results\TempPSD\'];
        else
            temppath = [outputfolder 'Results/TempPSD/'];
        end

        if exist(temppath) == 0 %#ok<*EXIST>
            % No PSD results, so exit analysis
            return
        end

        %Extract file names from folder
        files = dir(temppath);
        fileNames = {};
        for j=1:length(files)
            fileNames{j} = files(j).name;
        end

        %Read contents of each csv file
        contents = {};
        disp('Reading saved PSD results files from "TempPSD" folder...');
        for j = 1:length(fileNames);
            if strcmp(fileNames{j},'.') == 0 && strcmp(fileNames{j},'..') == 0
                contents{j} = read_mixed_csv([temppath fileNames{j}],',');
            else
                contents{j} = {};
            end
        end
        if length(fileNames)<3
            return
        end
        
        %%  Combine PSD results into single file
        for j = 3:length(contents);
            if j == 3
                finalResults = contents{j};
            else
                try
                    finalResults = [finalResults contents{j}];
                catch
                   msgbox('Failed to append PSD results into a single file.  Please remove your PSD results from the "Results\TempPSD" results folder before trying another analysis.') 
                end
            end
        end
        
        %%  Print combined results to file
        if ispc
            temppath2 = [outputfolder 'Results\'];
        else
            temppath2 = [outputfolder 'Results/'];
        end
        
        fid = fopen([temppath2 'Results(PSD)_' timeStamp '.csv'],'w');
        for j = 1:size(finalResults,1);
            for w = 1:size(finalResults,2);
                fprintf(fid,'%s',finalResults{j,w});
                fprintf(fid,'%s',',');
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        disp('...PSD results are now saved to a single csv file in the "Results" folder');
        
        %% Delete temp PSD files
        for j = 1:length(fileNames);
            if strcmp(fileNames{j},'.') == 0 && strcmp(fileNames{j},'..') == 0
                delete([temppath fileNames{j}]);
            end
        end
    end
        
    function lineArray = read_mixed_csv(fileName,delimiter)
      fid = fopen(fileName,'r');   % Open the file
      lineArray = cell(100,1);     % Preallocate a cell array (ideally slightly
                                   %   larger than is needed)
      lineIndex = 1;               % Index of cell to place the next line in
      nextLine = fgetl(fid);       % Read the first line from the file
      while ~isequal(nextLine,-1)         % Loop while not at the end of the file
        lineArray{lineIndex} = nextLine;  % Add the line to the cell array
        lineIndex = lineIndex+1;          % Increment the line index
        nextLine = fgetl(fid);            % Read the next line from the file
      end
      fclose(fid);                 % Close the file
      lineArray = lineArray(1:lineIndex-1);  % Remove empty cells, if needed
      for iLine = 1:lineIndex-1              % Loop over lines
        lineData = textscan(lineArray{iLine},'%s',...  % Read strings
                            'Delimiter',delimiter);
        lineData = lineData{1};              % Remove cell encapsulation
        if strcmp(lineArray{iLine}(end),delimiter)  % Account for when the line
          lineData{end+1} = '';                     %   ends with a delimiter
        end
        lineArray(iLine,1:numel(lineData)) = lineData;  % Overwrite line data
      end
    end
    function y = rmsX(x, dim)
        if nargin==1
          y = sqrt(mean(x .* conj(x)));
        else
          y = sqrt(mean(x .* conj(x), dim));
        end
    end
    function printResults(type)
        if writeOutputToFile == 1
            if ispc
                pathtemp = [outputfolder 'Results\Results(' type ')_' timeStamp '.csv'];
            else
                pathtemp = [outputfolder 'Results/Results(' type ')_' timeStamp '.csv'];
            end
            
            if exist(pathtemp) == 0;
                resultsFile = fopen(pathtemp,'w');
                fprintf(resultsFile,'%s\n',resultsHeader);
                fprintf(resultsFile,'%s\n',resultsSubHeader);
            else
                resultsFile = fopen(pathtemp,'a');
            end
            fprintf(resultsFile,'%s\n',results);
            fclose(resultsFile);
        end
    end
    function [f,Xmean,Xmedian,X5,X95] = PM_Analysis_fft(channel,window,olap,nfft,Fs)
        % olap = integer between 0 and 100
        
        %% Use global variables to save memory
        % trim data = trimmed waveform data
        switch channel
            case 'A'
                % Split data into windows segments
                data2ch1 = buffer(trimDataA,nfft,ceil(nfft*olap*1e-2),'nodelay'); % (allows overlap)
            case 'B'
                data2ch1 = buffer(trimDataB,nfft,ceil(nfft*olap*1e-2),'nodelay');
            case 'C'
                data2ch1 = buffer(trimDataC,nfft,ceil(nfft*olap*1e-2),'nodelay');
            case 'D'
                data2ch1 = buffer(trimDataD,nfft,ceil(nfft*olap*1e-2),'nodelay');
            case 'Aa'
                data2ch1 = buffer(trimDataAacc,nfft,ceil(nfft*olap*1e-2),'nodelay');
            case 'Ba'
                data2ch1 = buffer(trimDataBacc,nfft,ceil(nfft*olap*1e-2),'nodelay');
            case 'Ca'
                data2ch1 = buffer(trimDataCacc,nfft,ceil(nfft*olap*1e-2),'nodelay');
        end

        %% Define window type and constants
        switch windowS{1}
            case 'Hann';
                window = hann(windowL);
                windowScalingFactor = 0.54;
                noisePowerBandwidth = 1.36;
            case 'Hamming';
                window = hamming(windowL);
                windowScalingFactor = 0.5;
                noisePowerBandwidth = 1.5;
        end
            
        %% FFT
        %Apply scaling factor
        [~,n] = size(data2ch1);
        data2ch1 = data2ch1.*repmat(window,1,n)/windowScalingFactor;
        %fft
        fft_data_ch1 = abs(fft(data2ch1))./nfft;
        fft_amp_ch1=fft_data_ch1(1:nfft/2+1,:).^2;
        % Correction for the noise power bandwidwith
        fft_amp_ch1=2*fft_amp_ch1/noisePowerBandwidth;
        % Deriving the frequency units
        f = (Fs/2*linspace(0,1,nfft/2+1))';
        % Normalizing for 1 second
        wintime = nfft/Fs;      %time in seconds of window
        acc1 = fft_amp_ch1*wintime;
        
        %%  Calculate results
        % Mean
        Xmean = mean(acc1,2);    
        % Percentiles
        X5 = prctile(acc1,5,2);
        Xmedian = prctile(acc1,50,2);
        X95 = prctile(acc1,95,2);
        
        %%  Save window segments to file for append PSD option
        if exportWindowSegments == 1;
            if ispc()
                tempDir = [outputfolder 'WindowSegments\'];
            else
                tempDir = [outputfolder 'WindowSegments/'];
            end
            warning('off','all')
                mkdir([outputfolder 'WindowSegments']);
            warning('on','all')
            save([tempDir channel '_' num2str(sampleNum) '.mat'],'acc1','f','sampleStart','sampleEnd');
        end
    end

    function setCurrentFigure(fig)
        
        fignumbers = cell2mat(get(findall(0,'type','figure'), 'Number'));
        if any(ismember(fignumbers, fig))
             % else, reselect figure
            set(0, 'currentfigure', fig);
        else    
            % if no figure created, make new figure object
            figure(fig);
            set(fig,'defaultAxesFontName', fontType); % set figure font type
            set(fig,'DefaultAxesFontSize',fontSize); % set default font size
            set(fig,'color','w'); % change backround color to white
        end
        setFigureSize;
    end

    function output = consistencyAnalysis(data,dBthreshold)
        temp = sum(data < -1*(10^(dBthreshold/20))); %count num of elements less than thresh
        temp = temp + sum(data > (10^(dBthreshold/20))); %count num of elements greater than thresh
        output = (temp/length(data))*100; % convert to percent
    end

    function plotCalibrationValues(figh,values,interpolated)
    setCurrentFigure(figh);
       % draw calib values
       plot(interpolated{4}(:,1),interpolated{4}(:,2));
       hold on;
       plot(interpolated{1}(:,1),interpolated{1}(:,2),'r');
       plot(interpolated{2}(:,1),interpolated{2}(:,2),'g');
       plot(interpolated{3}(:,1),interpolated{3}(:,2),'m');
       hold off;

       xlabel('Frequency (Hz)');
       ylabel('Reciever sensitivity (dB re V/\muPa or V/(nm/s))');
       legend('H','X','Y','Z');
    end
end