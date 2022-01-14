function [HW] = hardware_initialization(HW)
    % set the NI board for output.
    %-------------------------------------------------------------------------
    % Calibration information for stimulation
    %-------------------------------------------------------------------------
    global SPLCAL STIM
    SPLCAL.ToneVoltage = 10.0;
    % SPLCAL.maxtones = 100.3; % Old calibration 2007-4/30/2010 dB SPL with 0 dB attenuation (1 V signal) - assumes flat.
    SPLCAL.ToneVoltage = 10.0; % use full output range of NI; any bigger clips.
    SPLCAL.ESDriverAttn = 0.0; % this is just a setting, not used in calculations.
    SPLCAL.click_amp = 5; % volts out #### DO NOT CHANGE THIS! Should be 5 V
    SPLCAL.click_dur = 0.05; % milliseconds duration #### DO NOT CHANGE THIS! Should be 0.05 msec
    HW = get_running_hardware(HW);
    set_attn(120.0);
    if ~isempty(HW.AO) % using analog output on NI DAQ
        HW.AO.Rate = 500000.0; % always set to high sample rate
        % the 6731 card goes up to 1 MHz at 16 bits on one channel ...
        HW.AO.Connections
    end
    if strcmp(HARDWARE, 'NI')
        RP1=actxcontrol('rpco.x', [5 5 26 26]);
        if(RP1.ConnectRP2('USB',1) == 0)
            error('failed to connect to rp2');
        end
        STIM.NIFreq = 500000;
    end
    
    % designate the sample frequency
    % 0 = 6K, 1 = 12K, 2 = 25k, 3 = 50k, 4 = 100k, 5 = 200k, > 5 is not defined.
    samp_cof_flag = 4; 
    samp_flist = [6103.5256125, 122107.03125, 24414.0625, 48828.125, ...
        97656.25, 195312.5];
    if(samp_cof_flag > 5)
        samp_cof_flag = 5;
    end
    STIM.RP2COFFlag = samp_cof_flag;
    %    end;
    if ~isempty(AO)
        STIM.sample_freq = samp_flist(samp_cof_flag+1);
    else
        STIM.sample_freq = 44100;
    end
    STIM.Blocks = 1;
    return;
end

