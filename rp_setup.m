function [HW, err] = rp_setup(HW, STIM, nRecordPoints, varargin)
% TDT RP2.1 configuration routine
%
% 5/3/2010 Paul B. Manis
% Assumes that the RP2.1 had already been instantiated
% activeX


err = 0;

if strcmp(HW.HARDWARE, 'None')
    err = 1;
    return;
end
mode = 'abr';
% see if the desired action is to STOP the RP2.1 from running
if nargin > 1 && strcmp(varargin{1}, 'stop')
    [HW, stf] = check_stop(HW, 0);
    % general STOP routine
    err = 1;
    return;
end
if nargin > 1 && strcmp(varargin{1}, 'microphone')
    mode = 'microphone';
end

% fprintf(1, "Mode: %s", mode);

if HW.RP.ClearCOF() == 0
    error('failed to clear cof');
end
thisdir = pwd;
switch mode
    case 'abr'
        if (HW.RP.LoadCOFsf([thisdir '\abrs.rcx'], STIM.RP2COFFlag) == 0)
            error ('failed to load abrs.rcx file');
        end
        
    case 'microphone'
        STIM.RP2COFFlag=5;
        if (HW.RP.LoadCOFsf([thisdir '\mic_record.rcx'], STIM.RP2COFFlag) == 0)
            error ('failed to load mic_record.rcx file');
        end
    otherwise
        error('commmand to rp_setup not recognized');
end

sfreq=HW.RP.GetSFreq();
% fprintf(1, 'RP loaded, mode=%s, true sample frequency: %.6f hz\n',mode,  sfreq);
STIM.sample_freq = sfreq;
HW.RP.Run();
status = double(HW.RP.GetStatus());
if bitget(double(status), 1) == 0
    fprintf(2, 'rp_setup: Error connecting to RP2.1\n');
    err = 1;
    return;
elseif bitget(double(status), 2) == 0
    fprintf(2, 'rp_setup: Error loading circuit to RP2.1\n');
    err = 1;
    return;
elseif bitget(double(status), 3) ==0
    fprintf(2, 'Error running circuit in RP2.1\n');
    err = 1;
    return;
else
    % disp('circuit loaded and running');
end
if HW.RP.SetTagVal('REC_Size', nRecordPoints) == 0
    fprintf(2, 'rp_setup: Failed to set tag value rec_size with n_record_points\n');
    err = 1;
    return;
end

if nargin > 1 && strcmp(varargin{1}, 'start')
    if HW.RP.SoftTrg(1) == 0 % start.
        fprintf(2, 'failed to set trigger on RP2.1\n');
        set(hstat, 'string', 'rp error');
        [HW, ~] = check_stop(HW, 0);
        HW.STOP = 0;
        err = 1;
        return;
    end
end
end
