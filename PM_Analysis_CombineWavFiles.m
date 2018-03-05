function PM_Analysis_CombineWavFiles

    answer = questdlg('Make a single multi-channel wav file or batch proccess to create many multi-channel wav files?','Processing options','Single','Batch','Single');

    % load saved paths
    try
        load('pathCalib.mat');
    catch
       pathCalib = '';
       pathFile = '';
       pathTone = '';
       fileName = '';
       calibFile = '';
    end

    switch answer
        case 'Single'
            %% Get locations of single chan wav files
            loop = true;
            channels = 0; % How many channels for output file
            chans = {'X' 'Y' 'Z' 'H'};
            paths = {}; % cell array of input paths
            while loop == true;
                if channels == 3;
                    loop = false;
                end

                %% Add extra channel
                answer = 'Yes';
                if channels > 0
                    answer = questdlg(['Do you want to add a ' chans{channels + 1} ' channel?'],'Add Channel','Yes','No','Yes');
                end

                %% Get file path
                if strcmp('Yes',answer);
                    [file, pathTone] = uigetfile({'*.wav' 'wav file'},['Select single channel wav file for channel ' chans{channels + 1} ],pathTone);
                    channels = channels + 1;
                    paths{channels} = [pathTone '/' file];
                    save('pathCalib.mat','pathCalib','pathFile','pathTone','fileName','calibFile');
                else
                    break;
                end 
            end
            
            %% Output
            path = uigetdir('','Select folder for output');
            if ispc
                path = [path '\'];
            else
                path = [path '/'];
            end
            name = inputdlg('Enter name for wav file output');
            if isempty(name)
                name = 'combined';
            else
                name = name{1};
            end
            combineWavs(paths, channels, [path name]);
        case 'Batch'
            %% Batch analyse wav files
            
            dirpath = '';
            try
                load('multichanpath');
            catch
            end
            
            disp('----------------------------');
            disp('-Create multi channel files-');
            
            % get file path
            dirpath = uigetdir(dirpath,'Select folder containing single-channel wav files');
            
            % exit script if no path was selected
            if dirpath == 0;
                return;
            end
            
            if ispc
                dirpath = [dirpath '\'];
            else
                dirpath = [dirpath '/'];
            end
            
            save('multichanpath','dirpath');
            
            %% Load files in folder
            if ispc
                files = dir([dirpath '*.wav']);
            else
                files = dir([dirpath '*.wav']);
            end
            
            %% User to input parsing mask
            uiwait(PM_Analysis_CombineWavFiles_batch);
            
            mask = '';
            A = '';
            B = '';
            C = '';
            D = '';
            
            load('batchProcessing.mat')
            
            disp(['Input mask:' mask]);
            %% Sort files into cell array
            fileNames = {};
            for i=1:length(files);
                fileNames{i,1} = files(i).name;
                %Parse file name
                Imask = 1;
                Ifile = 1;
                
                % reset names
                chan = '';
                fileName = '';
                
                while Imask <= length(mask) && Ifile <= length(fileNames{i,1});
                    switch mask(Imask)
                        case 'x'
                            Imask = Imask + 1;
                            Ifile = Ifile + 1;
                        case 'X'
                            %Move to next delimiter
                            delim = mask(Imask + 1);
                            while Ifile <= length(fileNames{i,1});
                                if fileNames{i,1}(Ifile) == delim
                                    break
                                end
                                Ifile = Ifile + 1;
                            end
                            Imask = Imask + 1;
                        case 'n'
                            fileName = [fileName fileNames{i,1}(Ifile)];
                            Imask = Imask + 1;
                            Ifile = Ifile + 1;
                        case 'N'
                            %Move to next delimiter
                            delim = mask(Imask + 1);
                            while Ifile <= length(fileNames{i,1});
                                if fileNames{i,1}(Ifile) == delim
                                    break
                                end
                                fileName = [fileName fileNames{i,1}(Ifile)];
                                Ifile = Ifile + 1;
                            end
                            Imask = Imask + 1;
                        case 'c'
                            chan = [chan fileNames{i,1}(Ifile)];
                            Imask = Imask + 1;
                            Ifile = Ifile + 1;
                        case 'C'
                            %Move to next delimiter
                            delim = mask(Imask + 1);
                            while Ifile <= length(fileNames{i,1});
                                if fileNames{i,1}(Ifile) == delim
                                    break
                                end
                                chan = [chan fileNames{i,1}(Ifile)];
                                Ifile = Ifile + 1;
                            end
                            Imask = Imask + 1;
                        otherwise
                            %Skip this caharcter, is a delimiter
                            Imask = Imask + 1;
                            Ifile = Ifile + 1;
                    end
                end
                if Imask > length(mask)
                    %Parse was succesfull!
                    fileNames{i,2} = chan;
                    fileNames{i,3} = fileName;
                else
                    %Parse failed
                    fileNames{i,2} = '';
                    fileNames{i,3} = '';
                end
            end

            for i=1:size(fileNames,1)
                switch fileNames{i,2}
                    case A
                        if strcmp(A,'') == 0
                            fileNames{i,2} = 'X';
                        end
                    case B
                        if strcmp(B,'') == 0
                            fileNames{i,2} = 'Y';
                        end
                    case C
                        if strcmp(C,'') == 0
                            fileNames{i,2} = 'Z';
                        end
                    case D
                        if strcmp(D,'') == 0
                            fileNames{i,2} = 'H';
                        end
                    otherwise
                        fileNames{i,2} = '';
                end
            end
            
            if ispc
                path = [dirpath '\Combined\'];    
            else
                path = [dirpath '/Combined/'];
            end
            warning('off','all');
                mkdir(path);
            warning('on','all');
            %% Group files and merge into 
            filesLabels = unique(fileNames(:,3));  
              
            filesLabesls2 = {};
            count = 1;
            disp('Unique identifiers:');
            for j = 1:length(filesLabels)
                if strcmp(filesLabels{j},'') == 0
                    filesLabesls2{count} = filesLabels{j};
                    count = count + 1;
                    disp(filesLabels{j});
                end
            end
            
            filesLabels = filesLabesls2;
            
            disp('-----Processing files-----');
            
            if length(filesLabels) == 0;
                disp('No files identified which match your parsing criteria...');
                disp('...stopping script');
            end
            
            mergeCount = 0;
            for i = 1:length(filesLabels)
                disp(['File:' filesLabels{i} '...']);
                %Find indexes of matched file names
                xf = 0;
                yf = 0;
                zf = 0;
                hf = 0;
                for j=1:size(fileNames,1);
                    if strcmp(fileNames{j,2},'X') > 0 && strcmp(filesLabels(i),fileNames{j,3})
                        xf = j;
                        break;
                    end
                end
                if xf > 0
                    disp('...X channel detected');
                    for j=1:size(fileNames,1);
                        if strcmp(fileNames{j,2},'Y') > 0 && strcmp(filesLabels(i),fileNames{j,3})
                            yf = j;
                            break;
                        end
                    end
                    if yf > 0
                        disp('...Y channel detected');
                        for j=1:size(fileNames,1);
                            if strcmp(fileNames{j,2},'Z') > 0 && strcmp(filesLabels(i),fileNames{j,3})
                                zf = j;
                                break;
                            end
                        end
                        if zf > 0
                            disp('...Z channel detected');
                            for j=1:size(fileNames,1);
                                if strcmp(fileNames{j,2},'H') > 0 && strcmp(filesLabels(i),fileNames{j,3})
                                    hf = j;
                                    break;
                                end
                            end
                            if hf > 0
                                disp('...H channel detected');
                                %merge XYZH
                                if ispc
                                    combineWavs({[dirpath '\' fileNames{xf,1}];[dirpath '\' fileNames{yf,1}];[dirpath '\' fileNames{zf,1}];[dirpath '\' fileNames{hf,1}]},4,[path filesLabels{i}]);
                                else
                                    combineWavs({[dirpath '/' fileNames{xf,1}];[dirpath '/' fileNames{yf,1}];[dirpath '/' fileNames{zf,1}];[dirpath '/' fileNames{hf,1}]},4,[path filesLabels{i}]);
                                end
                                mergeCount = mergeCount + 1;
                            else
                                %merge XYZ
                                if ispc
                                    combineWavs({[dirpath '\' fileNames{xf,1}];[dirpath '\' fileNames{yf,1}];[dirpath '\' fileNames{zf,1}]},3,[path filesLabels{i}]);
                                else
                                    combineWavs({[dirpath '/' fileNames{xf,1}];[dirpath '/' fileNames{yf,1}];[dirpath '/' fileNames{zf,1}]},3,[path filesLabels{i}]);
                                end
                                mergeCount = mergeCount + 1;
                            end
                        else
                            %merge XY
                            if ispc
                                combineWavs({[dirpath '\' fileNames{xf,1}];[dirpath '\' fileNames{yf,1}]},2,[path filesLabels{i}]);
                            else
                                combineWavs({[dirpath '/' fileNames{xf,1}];[dirpath '/' fileNames{yf,1}]},2,[path filesLabels{i}]);
                            end
                            mergeCount = mergeCount + 1;
                        end
                    else
                        disp('...No Y channel detected')
                    end
                else    
                    disp('...No X channel detected')
                end
            end
            disp('-----Operation complete-----');
            disp([num2str(mergeCount) ' files created from ' num2str(length(filesLabels)) ' unique identifiers']);
            beep;
            disp(['Processed files saved at this location:' path]);
            
            disp('All created wav files are 16-bit');
            clear all;
            clear global
            warning('off','all');
                delete 'TempCSVDataA.dat' 'TempCSVDataB.dat' 'TempCSVDataC.dat' 'TempCSVDataD.dat';
            warning('on','all');
    end
    
    function combineWavs(paths,channels,name)
        disp(['...creating ' num2str(channels) ' channel wav file'])
        if channels > 1;
            % Create files to save wavdata from origional files to memmap variable
            fidWriteA = fopen('TempCSVDataA.dat','w');
            fidWriteB = fopen('TempCSVDataB.dat','w');
            fidWriteC = fopen('TempCSVDataC.dat','w');
            fidWriteD = fopen('TempCSVDataD.dat','w');

            % Get wav file size (needed to know how much hardd drive space we need to store data!)
            if verLessThan('matlab', '8.1');
                [temp,samplez] = wavfinfo(paths{1});
                samplez = str2double(samplez(29:length(samplez)-24));
                wavInfo(1).TotalSamples = samplez;
                [y,fs] = wavread(paths{1},1);
                wavInfo(1).NumChannels = size(y,2);
                wavInfo(1).SampleRate = fs;

                [temp,samplez] = wavfinfo(paths{2});
                samplez = str2double(samplez(29:length(samplez)-24));
                wavInfo(2).TotalSamples = samplez;
                [y,fs] = wavread(paths{1},1);
                wavInfo(2).NumChannels = size(y,2);
                wavInfo(2).SampleRate = fs;
                if channels > 2
                    [temp,samplez] = wavfinfo(paths{3});
                    samplez = str2double(samplez(29:length(samplez)-24));
                    wavInfo(3).TotalSamples = samplez;
                    [y,fs] = wavread(paths{1},1);
                    wavInfo(3).NumChannels = size(y,2);
                    wavInfo(3).SampleRate = fs;
                end
                if channels > 3
                    [temp,samplez] = wavfinfo(paths{4});
                    samplez = str2double(samplez(29:length(samplez)-24));
                    wavInfo(4).TotalSamples = samplez;
                    [y,fs] = wavread(paths{1},1);
                    wavInfo(4).NumChannels = size(y,2);
                    wavInfo(4).SampleRate = fs;
                end
            else
                wavInfo(1) = audioinfo(paths{1});
                wavInfo(2) = audioinfo(paths{2});
                if channels > 2
                    wavInfo(3) = audioinfo(paths{3});
                end
                if channels > 3
                    wavInfo(4) = audioinfo(paths{4});
                end
            end

            %% Check that wav files all have 1 channel
             check = 0;
             for i=1:length(wavInfo)
                 check = max(check,wavInfo(i).NumChannels);
             end

             if (check > 1)
                 disp('At least one of the selected wav files has more than 1 channel!!!');
                 disp('You may only combine single channel wav files.')
                 return;
             end
            %% Save smallest wavfile length
            check = wavInfo(1).TotalSamples;
            for i=1:length(wavInfo)
                check = min(check,wavInfo(i).TotalSamples);
            end
            samples = check;

            %% Check that all wavs have same samplerate
            % save samllest wavfile length
            checkmin = wavInfo(1).SampleRate;
            checkmax = 0;
            for i=1:length(wavInfo)
                checkmin = min(checkmin,wavInfo(i).SampleRate);
                checkmax = max(checkmax,wavInfo(i).SampleRate);
            end
            if ~(checkmin == checkmax)
                disp('Selected files do not have the same sample rate');
                return;
            end
            
            disp('...copying wav files to temporary containers...')
            %% Start copying data to temp files
            samplesIndex = 1;
            bufferSize = 400000; %reduce this if your having "out of memory errors"

            while (samplesIndex<samples)
                if bufferSize + samplesIndex > samples
                    linesToWrite = int32(samples - samplesIndex);
                else
                    linesToWrite = int32(bufferSize);
                end

                % Load portion of data from input wav files and convert to
                % 16int
                [tempA,SampleRate] = audioread(paths{1}, double([samplesIndex (samplesIndex + linesToWrite)]),'native');
                tempA = int16(tempA);
                [tempB,SampleRate] = audioread(paths{2}, double([samplesIndex (samplesIndex + linesToWrite)]),'native');
                tempB = int16(tempB);
                if channels > 2
                    [tempC,SampleRate] = audioread(paths{3}, double([samplesIndex (samplesIndex + linesToWrite)]),'native');
                    tempC = int16(tempC);
                end
                if channels > 3
                    [tempD,SampleRate] = audioread(paths{4}, double([samplesIndex (samplesIndex + linesToWrite)]),'native');
                    tempD = int16(tempD);
                end

                % write to temp files in harddrive
                fwrite(fidWriteA,tempA(:),'int16');
                fwrite(fidWriteB,tempB(:),'int16');
                if channels > 2
                    fwrite(fidWriteC,tempC(:),'int16');
                end
                if channels > 3
                    fwrite(fidWriteD,tempD(:),'int16');
                end
                samplesIndex = int32(samplesIndex + linesToWrite);
            end
            
            fclose(fidWriteA);
            fclose(fidWriteB);
            fclose(fidWriteC);
            fclose(fidWriteD); 

            %% write data to final wav file
            disp('...writing final multi-chan wav file...')
            if ispc
                path = [path '\'];
            else
                path = [path '/'];
            end

            filename = [name '_' num2str(channels) 'channels.wav'];
            loopS = true;
            %% Write to wav file (recursevly limit file size if memory error occurs)
            switch channels
                case 1
                    disp('Not enough channels in file')
                    return
                case 2
                    PM_Analysis_wavwriteMemmap({'TempCSVDataA.dat' 'TempCSVDataB.dat'},SampleRate,filename)
                case 3
                    PM_Analysis_wavwriteMemmap({'TempCSVDataA.dat' 'TempCSVDataB.dat' 'TempCSVDataC.dat'},SampleRate,filename)
                case 4
                    PM_Analysis_wavwriteMemmap({'TempCSVDataA.dat' 'TempCSVDataB.dat' 'TempCSVDataC.dat' 'TempCSVDataD.dat'},SampleRate,filename)
            end
                    
            disp([num2str(channels) ' channel wav file saved as ' filename]);
            beep;
        end
    end
end