function Calibration_Test
    SampleRate = 96000
    maxCalibrationVal = SampleRate/2
    minCalibrationVal = 1

    %Create array with calibration values and frequencies
    % I chose the arbitrary value of -80 V/uPa for the calibration array
    interpConFD = 1:20:48000;
    interpConFD(2,:) = -80;
  
    %import waveform
    RawDataD.data = wavread('Impulse_Test.wav');
    
    cConstants = interpConFD';
    %[f0s,Spect]=fftFunc(RawDataD.data,1/SampleRate);
    [f0s,Spect]=fftFunc(RawDataD.data,1/SampleRate);

    % Equation for converting Omni channel
    % into SPL for testing
    % SPL = 20log(V)+176
    % 20*log10(rms(RawDataD.Data(0.2e5:1.3e5))) + 176
    % Use this equation to proof your results to make sure the
    % calibration equations are correct

    %%  Interpolate calibration values and convert to linear units
    %Interpolate reciever sensitivity values
    %from calibration data

    %Interpolate missing values in calibration data
    interpCalib = spline(cConstants(:,1),20*log10(cConstants(:,2)),f0s);
    interpCalib = 10.^(interpCalib./20); 

    % Write flat line on calibration values outside
    % calibration range
    interpCalib(1:(find(f0s>=cConstants(1,1),1,'first') )) = cConstants(1,2);
    interpCalib((find(f0s>=cConstants(end,1),1,'first') ):length(f0s)) = cConstants(end,2);

%% Convert from V to Pa(or m/s)
    %Convert raw data V into uPa (or m/s) for each
    %frequency
    Spect = Spect./(interpCalib.');

%% IFFT (convert the signal back into time domain) and save to hard drive
    % Update new sample rate for calibrated data
    Fs = 2*f0s(end);
    SampleRate = Fs;

    % Create padding for data
    NFFT =(2.*length( Spect))-1; %NFFT points(Next power of 2 from length)

    % Convert signal back into time domain (xt)
    xt = (1/sqrt(2))*NFFT*ifft(Spect,NFFT,'symmetric');
    % Filter frequencies outside opf
    % calibration range
    % Write Calibrated Data to harddrive and
    % reassign memmap variable


    %%  Plot results
    
    % X axis limits
    a = 3e4;
    b = 3.5e4;

    figure(1);
    subplot(2,1,1)
    plot(real(RawDataD.data));
    grid on
    %xlim([a,b])
    subplot(2,1,2)
    plot(real(xt),'r');
    %xlim([a,b])
    grid on

    figure(2);
    plot(real(RawDataD.data),'LineWidth',2);
    hold on
    % plot normalized wave functions against eachother
    plot((real(-xt)/max(real(xt)))*max(real(RawDataD.data)),'r');
    hold off
    %xlim([a,b])
    grid on
    
                               
    function [freq, P]=fftFunc(p,dt)
        % calculates FFT of given signal
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
     % Reenable all warnings
        warning('on','all');
        freq = (0.5/dt)*linspace(0,1,NFFT/2+1); % frequency vector

        % NFFT = 2^nextpow2(nsamp); % Next power of 2 from length of y
        % Y= fft(p,NFFT)/nsamp;
        % freq = (1/dt)/2*linspace(0,1,NFFT/2);
        % 
        % P=2*abs(Y(1:NFFT/2)); 
    end

    function PM_Calibration_iFFT
        %% FFT 
        [ f0s Spect]=fftFunc(input_signal,dt);
        phase_data= exp(i*angle(Spect));

        %% IFFT
        Fs =2*f0s(end);
        NFFT    =(2.*length( Spect))-1; %#FFT points(Next power of 2 from length)
        xt = (1/sqrt(2))*NFFT*ifft( Spect,NFFT,'symmetric');
        T = [ 0 : NFFT-1] / Fs;
    end
end