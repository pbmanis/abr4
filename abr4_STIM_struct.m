classdef abr4_STIM_struct
    % stim properties for abr4
    % holds all of the hardware structure information for abr4
    
    properties
        % Hardware structure
        Info  % text string identifying this structure
        RP2COFFlag  % RP2 flag (does this belong in HW?)
        period  % stimulus period (time between stimuli)
        delay  % stimulus delay in msec from trigger
        ipi  % interpulse interval for repeated stimuli
        risefall  % rise fall time for tone pips, in msec
        NIFreq  % NI output frequency.
        sample_freq  % RP2.1 input sample frequency
        rate  % 1./sample_freq
        Blocks  % 
        tone_delay  % delay for tone pip
        tone_delay_mapping  % delay for tone pip in mapping mode
        click_delay  % delay for click
        click_dur  % duration of each click
        click_amp  % amplitude of click pulse
        StimPerSweep  % number of click stimuli or tone pips per sweep
        Alternate  % flag for alternation of click sign
        NSweeps  % number of sweeps 
        InterSweepInterval  % interval between sweeps
        avg_dur  % 
        avg_reject
        rms_reject
        amp_gain
        SPLMin
        SPLMax
        SPLStep
        spls  % spls in the run
        freqs  % frequencies in the run
        FreqList  % string describing the frequency list
        maxtones  
        maxclick
        wave  % acoustic stimulus waveform
        clock  % clock for NI waveform.
        Monitor
        ACQPars_np  % returned by acquire 4... 
        ACQPars_nn
        ACQPars_tb
        ACQPars_nall
    end
    
    methods
        function obj = initialize(obj)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Info = 'ABR4 StimFile';
            obj.RP2COFFlag = [];
            obj.sample_freq = 10000.0;
            obj.rate = 1e-4;
            obj.Blocks = 1;
            obj.period = 1.0;
            obj.delay = 1.0;
            obj.ipi = 10.;
            obj.risefall = 0.1;
            obj.NIFreq = 500000;
            obj.tone_delay = 0.1;
            obj.tone_delay_mapping = 10.;
            obj.click_delay = 1.0;
            obj.click_dur = 0.1;
            obj.click_amp = 5.0;
            obj.StimPerSweep = 1;
            obj.Alternate = 0;
            obj.hsr = 1.0;
            obj.hisi = 40.0;
            obj.hrep = 1.0;
            obj.NSweeps = 100.;
            obj.InterSweepInterval = 1.5;
            obj.avg_dur = 1.0;
            obj.avg_reject = 0;
            obj.rms_reject = 0;
            obj.amp_gain = 10000;
            obj.spls = [0, 50];
            obj.SPLMin = 0;
            obj.SPLMax = 90;
            obj.SPLStep = 10;
            obj.freqs = [];
            obj.FreqList = '';
            obj.wave = [];
            obj.clock = 0.;
            obj.Monitor = [];
            obj.ACQPars_np = 0;
            obj.ACQPars_nn = 0;
            obj.ACQPars_tb = [];
            obj.ACQPars_nall = 0;
        end
    end
end

