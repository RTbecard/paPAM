function combinePSDs(outputfolder,FileName,timeStamp)
    f = [];
    acc1 = [];
    sampleStart = [];
    sampleEnd = [];
    timeRanges = '';
    %% Load files from folder
    if ispc
        folder = [outputfolder 'WindowSegments\'];
    else
        folder = [outputfolder 'WindowSegments/'];
    end

    files = dir(folder);
    fileNames = {};
    for i = 1:length(files)
        if strcmp(files(i).name,'.')   == 0 && strcmp(files(i).name,'..') == 0
            fileNames = {fileNames{:} files(i).name};
        end
    end

    channelNames = {'A' 'B' 'C' 'D' 'Aa' 'Ba' 'Ca'};

    loops = 1;
    
    resultsWindow = {};
    
    %% Loop channel names
    for q = 1:length(channelNames)
        appendedWindows = [];
        I = strfind(fileNames,[channelNames{q} '_']);
        I = find(cellfun(@isempty,I)==0);
        %% Loop files within channel name and append data
        for i=1:length(I)
            load([folder fileNames{I(i)}]);
            timeRanges = [timeRanges num2str(sampleStart) '-' num2str(sampleEnd) ' '];
            appendedWindows = [appendedWindows acc1];
        end
        %% Save results in files
        if isempty(I) == 0;
            Xmean = mean(acc1,2);
            X5 =  prctile(acc1,5,2);
            X50 = prctile(acc1,50,2);
            X95 = prctile(acc1,95,2);

            resultsWindow{1,(loops-1)*5 + 1} = channelNames{q};
            resultsWindow{1,(loops-1)*5 + 2} = ['n = ' num2str(size(appendedWindows,2))];
            resultsWindow{1,(loops-1)*5 + 3} = timeRanges;
            resultsWindow{2,(loops-1)*5 + 1} = 'Hz';
            resultsWindow{2,(loops-1)*5 + 2} = 'Mean';
            resultsWindow{2,(loops-1)*5 + 3} = '5th';
            resultsWindow{2,(loops-1)*5 + 4} = '50th';
            resultsWindow{2,(loops-1)*5 + 5} = '95th';
            for w = 1:length(f);
                resultsWindow{w+2,(loops-1)*5 + 1} = num2str(f(w));
            end
            for w = 1:length(f);
                resultsWindow{w+2,(loops-1)*5 + 2} = num2str(10*log10(Xmean(w)));
                resultsWindow{w+2,(loops-1)*5 + 3} = num2str(10*log10(X5(w)));
                resultsWindow{w+2,(loops-1)*5 + 4} = num2str(10*log10(X50(w)));
                resultsWindow{w+2,(loops-1)*5 + 5} = num2str(10*log10(X95(w)));
            end
            loops = loops + 1;
        end
    end

    %% Print results to file
    if ispc
        tempdir = [outputfolder 'Results\' 'appendedPSD_' timeStamp '_' FileName '.csv'];
    else
        tempdir = [outputfolder 'Results/' 'appendedPSD_' timeStamp '_' FileName '.csv'];
    end
    fw = fopen(tempdir, 'w') ;
    for q = 1:size(resultsWindow,1)
        fprintf(fw, '%s,', resultsWindow{q,1:end-1});
        fprintf(fw, '%s\n', resultsWindow{q,end}) ;
    end
    fclose(fw);
    
    if ispc
        delete([outputfolder 'WindowSegments\*.mat']);
    else
        delete([outputfolder 'WindowSegments/*.mat']);
    end