function PM_Analysis_wavwriteMemmap(y,Fs,wavefile)
  
    % Parse inputs:
    error(nargchk(2,4,nargin));
    
    %Default to 16biut files
    nbits = 16;
 
    % get num of channels
    channels = length(y);

    %Get Length of Y
    for i  = 1:length(y)
       y{i} = memmapfile(y{i},'Writable',false,'Format','int16');
    end
    samples = length(y{1}.Data);

    % Determine number of bytes in chunks
    % (not including pad bytes, if needed):
    % ----------------------------------
    %  'RIFF'           4 bytes
    %  size             4 bytes
    %  'WAVE'           4 bytes
    %  'fmt '           4 bytes
    %  size             4 bytes
    % <wave-format>     14 bytes
    % <format_specific> 2 bytes (PCM)
    %  'data'           4 bytes
    %  size             4 bytes
    % <wave-data>       N bytes
    % ----------------------------------
    bytes_per_sample = ceil(nbits/8);
    total_samples    = samples * channels;
    total_bytes      = total_samples * bytes_per_sample;

    riff_cksize = 36+total_bytes;   % Don't include 'RIFF' or its size field
    fmt_cksize  = 16;               % Don't include 'fmt ' or its size field
    data_cksize = total_bytes;      % Don't include 'data' or its size field

    % Determine pad bytes:
    data_pad    = rem(data_cksize,2);
    riff_cksize = riff_cksize + data_pad; % + fmt_pad, always 0

    % Open file for output:
    fid = OpenWaveWrite(wavefile);

    % file is now open, wrap the rest of the calls
    % in a try catch so we can close the file if there is a failure
    try
        % Prepare basic chunk structure fields:
        ck=[]; ck.fid=fid; ck.filename = wavefile;

        % Write RIFF chunk:
        ck.ID   = 'RIFF';
        ck.Size = riff_cksize;
        write_ckinfo(ck);

        % Write WAVE subchunk:
        ck.ID   = 'WAVE';
        ck.Size = [];  % Indicate a subchunk (no chunk size)
        write_ckinfo(ck);

        % Write <fmt-ck>:
        ck.ID   = 'fmt ';
        ck.Size = fmt_cksize;
        write_ckinfo(ck);

        % Write <wave-format>:
        fmt.filename        = wavefile;
        if nbits == 32,
            fmt.wFormatTag  = 3;            % Data encoding format (1=PCM, 3=Type 3 32-bit)
        else
            fmt.wFormatTag  = 1;
        end
        fmt.nChannels       = channels;     % Number of channels
        fmt.nSamplesPerSec  = Fs;           % Samples per second
        fmt.nAvgBytesPerSec = channels*bytes_per_sample*Fs; % Avg transfer rate
        fmt.nBlockAlign     = channels*bytes_per_sample;    % Block alignment
        fmt.nBitsPerSample  = nbits;        % standard <PCM-format-specific> info
        write_wavefmt(fid,fmt);

        % Write <data-ck>:
        ck.ID   = 'data';
        ck.Size = data_cksize;
        write_ckinfo(ck);

        % Write <wave-data>, and its pad byte if needed:
        write_wavedat(fid,fmt);

        % Close file:
        fclose(fid);
    catch exception
        fclose(fid);
        throw(exception);
    end
    function write_wavedat(fid,fmt)
        if fmt.wFormatTag==1 || fmt.wFormatTag==3,
           % PCM Format

           %Data must be in 16bit units at this point! 
           dtype='int16'; % signed 16-bit

           %Divide data into smaller packets for writing
           packetSize = 30*(5e5); %n*5e5 = n Mb of space required
           packets = ceil(samples/packetSize);
           
           % Write data to file!
           for j=1:packets
               % Copy data to matrix
               temp = [];
               for c=1:channels
                   if j == packets
                       temp = [temp y{c}.Data(((j-1)*packetSize)+1:end)];
                   else
                       temp = [temp y{c}.Data(((j-1)*packetSize)+1:j*packetSize)];
                   end
               end
               % Write to file
               fwrite(fid, temp', dtype);
               disp(['...' num2str(floor(100*((channels-1)*packets + j)/(packets*channels))) '% done writing file...']);
           end

           % Determine # bytes/sample - format requires rounding
           %  to next integer number of bytes:
           BytesPerSample = ceil(fmt.nBitsPerSample/8);
           %nBitsPerSampe = 16 forced
           
           % Determine if a pad-byte must be appended to data chunk:
           if rem(total_samples*BytesPerSample, 2) ~= 0,
              fwrite(fid,0,'uint8');
           end

        else
          % Unknown wave-format for data.
          error(message('MATLAB:audiovideo:wavewrite:unsupportedDataFormat'));
        end

        return
    end
end

function [fid] = OpenWaveWrite(wavefile)
        % OpenWaveWrite
    %   Open WAV file for writing.
    %   If filename does not contain an extension, add ".wav"

    if ~ischar(wavefile),
       error(message('MATLAB:audiovideo:wavewrite:InvalidFileNameType')); 
    end
    if isempty(findstr(wavefile,'.')),
      wavefile=[wavefile '.wav'];
    end
    % Open file, little-endian:
    [fid,err] = fopen(wavefile,'wb','l');
    if (fid == -1)
        error(message('MATLAB:audiovideo:wavewrite:unableToOpenFile', err));
    end
    return
end

% ------------------------------------------------------------------------
function write_ckinfo(ck)
    % WRITE_CKINFO: Writes next RIFF chunk, but not the chunk data.
    %   Assumes the following fields in ck:
    %         .fid   File ID to an open file
    %         .ID    4-character string chunk identifier
    %         .Size  Size of chunk (empty if subchunk)
    %
    %
    %   Expects an open FID pointing to first byte of chunk header,
    %   and a chunk structure.
    %   ck.fid, ck.ID, ck.Size, ck.Data

    if (fwrite(ck.fid, ck.ID, 'char') ~= 4),
       error(message('MATLAB:audiovideo:wavewrite:failedChunkInfoWrite', num2str( ck.ID ), ck.filename));
    end

    if ~isempty(ck.Size),
      % Write chunk size:
      if (fwrite(ck.fid, ck.Size, 'uint32') ~= 1),
          error(message('MATLAB:audiovideo:wavewrite:failedChunkInfoWrite', num2str( ck.ID ), ck.filename));
      end
    end

    return
end

% ------------------------------------------------------------------------
function write_wavefmt(fid, fmt)
    % WRITE_WAVEFMT: Write WAVE format chunk.
    %   Assumes fid points to the wave-format subchunk.
    %   Requires chunk structure to be passed, indicating
    %   the length of the chunk.

    % Create <wave-format> data:
    if (fwrite(fid, fmt.wFormatTag,      'uint16') ~= 1) || ...
       (fwrite(fid, fmt.nChannels,       'uint16') ~= 1) || ...
       (fwrite(fid, fmt.nSamplesPerSec,  'uint32' ) ~= 1) || ...
       (fwrite(fid, fmt.nAvgBytesPerSec, 'uint32' ) ~= 1) || ...
       (fwrite(fid, fmt.nBlockAlign,     'uint16') ~= 1),
       error(message('MATLAB:audiovideo:wavewrite:failedWaveFmtWrite', fmt.filename));
    end

    % Write format-specific info:
    if fmt.wFormatTag==1 || fmt.wFormatTag==3,
      % Write standard <PCM-format-specific> info:
      if (fwrite(fid, fmt.nBitsPerSample, 'uint16') ~= 1),
         error(message('MATLAB:audiovideo:wavewrite:failedWaveFmtWrite', fmt.filename));
      end

    else
      error(message('MATLAB:audiovideo:wavewrite:unknownDataFormat'));
    end

    return
end
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------


% end of wavwrite.m
