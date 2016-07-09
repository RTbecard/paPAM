function PM_Analysis_CreateCalibTone()
    path = uigetdir('','Export clibration tone to folder');
 
    if length(path) < 2
        return
    end
 
    if ispc
        path = [path '\'];
    else
        path = [path '/'];
    end
    duration = inputdlg({'Duration of tone (seconds):' 'Frequency (Hz)' 'Amplitude (dBFS):'},'Tone Settings',1,{'10' '200' '-6'});
    
    if isempty(duration)
        return
    end
    
    fs = 44100;
    seconds = str2double(duration{1});
    freq = str2double(duration{2});
    amplitude = str2double(duration{3});
    amplitude = 10^(amplitude/20);
    ending = 2*pi*freq;
    
    y = (sin(0:ending/fs:ending*seconds)) * amplitude;
    
    dBFS = num2str(20*log10(amplitude));
    freqS = num2str(freq);
    
    name = ['calibrationTone_' dBFS 'dBFS_' freqS 'Hz.wav'];
    
    if verLessThan('matlab', '8.1');
        wavwrite(y,fs,[path name]);
    else
        audiowrite([path name],y,fs);
    end
    
    beep;
    msgbox(['Calibration tone has been saved to: ' path 'as ' name]);
end

