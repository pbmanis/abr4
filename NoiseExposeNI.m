function NoiseExposeNI(varargin)
% NOISEEXPOSE M-file for NoiseExpose.fig
%      NOISEEXPOSE, by itself, creates text new NOISEEXPOSE or raises the existing
%      singleton*.
%
%      H = NOISEEXPOSE returns the handle to text new NOISEEXPOSE or the handle to
%      the existing singleton*.
%
%      NOISEEXPOSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NOISEEXPOSE.M with the given input arguments.
%

global QUIT RP PA5 AO UseAttn MicIn
global Freq Band_Width AO_freq Attn Duration
global PlotHandles Hspec Hwave
persistent swave fs continueFlag remainingTime elapsed totalDuration total_elapsed
AO_freq = 500000; % (note this is default for noise_gen.m, and may need to be passed in the future)
UseAttn = 'right'; % set to 'left' to force use of left attenuator.
persistent topDBSPL

topDBSPL = 122.8; % dB SPL, Crown Channel 2 at 15 (green arrow), 0db Attn

if(nargin == 0) % initialize
    hf = findobj('Tag', 'NoiseExposeFig');
    if(~isempty(hf))
        fprintf(1, 'Only 1 Instance of NoiseExposeNI is allowed!\n');
        return;
    end;
    RP = [];
    PA5 = [];
    Freq = 12000;
    Band_Width = 4000;
    Attn = 120; % dB
    Duration = 2.0; % hours
    continueFlag = 0;
    remainingTime = 0;
    elapsed = 0;
    totalDuration = 0;
    total_elapsed = 0;
    open([pwd '\TFRNoiseExposeNIFig.fig']); % get the figure window
    figurehandle = gcf;
    PlotHandles.spectrum = findobj('tag', 'NoiseExposeGraph');
    PlotHandles.waveform = findobj('tag', 'waveform');    
    % inialize the GUI
    % Set the data from the GUI
    hcf=findobj('Tag', 'CF');
    if(~isempty(hcf) && ishandle(hcf))
        set(hcf, 'String', num2str(Freq));
    end;
    hbw=findobj('Tag', 'BW');
    if(~isempty(hbw) && ishandle(hbw))
        set(hbw, 'String', num2str(Band_Width));
    end;
    hat=findobj('Tag', 'AT');
    if(~isempty(hat) && ishandle(hat))
        set(hat, 'String', num2str(Attn));
    end;
    hed=findobj('Tag', 'ED');
    if(~isempty(hed) && ishandle(hed))
        set(hed, 'String', sprintf('%7.2f', Duration));
    end;
    hatn = findobj('Tag', 'AttnID');
    isempty(hatn)
    ishandle(hatn)
    if(~isempty(hatn) && ishandle(hatn))
        set(hatn, 'String', UseAttn);
    end;
    % access the NI Card
    daqreset;
%     get_running_hardware();
    setup_NI()
%    AO = analogoutput('nidaq', 'Dev1'); % for Matlab2006a
    %    AO = analogoutput('nidaq', 1); % forearlier versionts
%     set(AO, 'samplerate', AO_freq); % always set to high sample rate
    AO.Rate = AO_freq;
    % the 6731 card goes up to 1 MHz at 16 bits on one channel ...
    %set(AO, 'triggertype', 'immediate');        % compute the noise
%     addTriggerConnection(AO, 'External', 'Dev1/PFI0', 'StartTrigger');
%     Commented out above line on 3/9/17
%    AO.Connections
    return; % just jump back
end

% process button calls.

switch(varargin{1})
    case 'quit'
        if(ishandle(AO) && isvalid(AO))
            stop(AO);
            delete(AO);
        end;
        delete(gcf);
        return;
    case 'SetCF' % set the center frequency
        hcf=findobj('Tag', 'CF');
        if(~isempty(hcf) && ishandle(hcf))
            Freq = str2double(get(hcf, 'String'));
        end;
    case 'SetBW'
        hbw=findobj('Tag', 'BW');
        if(~isempty(hbw) && ishandle(hbw))
            Band_Width = str2double(get(hbw, 'String'));
        end;
    case 'SetLevel'
        hat=findobj('Tag', 'AT');
        if(~isempty(hat) && ishandle(hat))
            Attn = str2double(get(hat, 'String'));
        end;
    case 'SetDuration'
        hed=findobj('Tag', 'ED');
        if(~isempty(hed) && ishandle(hed))
            Duration = str2double(get(hed, 'String'));
        end;
        fprintf(2, 'New Duration: %7.2f hours\n', Duration);
        
    case {'preview', 'noisecalc'}
        [swave, fs] = noise_gen('bandpass', Freq-Band_Width/2.0, Freq+Band_Width/2.0, 'butter', 4,1);
       
        Hs = spectrum.yulear(1024);
        Hpsd = psd(Hs, swave,'Fs', 1.0/fs);
        plot(PlotHandles.spectrum, Hpsd.Frequencies(2:end)/1000., ...
                Hpsd.Data(2:end));
        set(PlotHandles.spectrum, 'XScale', 'log')
        set(PlotHandles.spectrum, 'YScale', 'log')
       set(PlotHandles.spectrum, 'XLim', [1. 100.]);
        %set(PlotHandles.spectrum, 'YLim', [-140 -20]);
         set(PlotHandles.spectrum, 'XTick', [1, 2, 5, 10, 20, 50, 100]);
         set(PlotHandles.spectrum, 'XTickLabel', {'1.0', '2.0', '5.0', '10.0', '20.0', '50.0', '100.0'});
        xlabel('kHz');

        tb = [0:fs:fs*(length(swave)-1)];
        plot(PlotHandles.waveform, tb, swave);
        
        
    case 'checklevel' % read the microphone and check the level/output
        load(sprintf('microphone_7016.cal'), '-mat'); % 1/2 " MICROPHONE
        [swave, fs] = noise_gen('bandpass', Freq-Band_Width/2.0, Freq+Band_Width/2.0, 'butter', 4,1);
        [ch1, ch2, AIsamp, err] = checklevel(swave, fs);
        
        Hs = spectrum.yulear(1024);
        Hpsd = psd(Hs, ch2, 'Fs', AIsamp);
        plot(PlotHandles.spectrum, Hpsd.Frequencies(2:end)/1000., ...
            Hpsd.Data(2:end));
        set(PlotHandles.spectrum, 'XScale', 'log')
        set(PlotHandles.spectrum, 'YScale', 'log')

        set(PlotHandles.spectrum, 'XLim', [1. 100.]);
        %set(PlotHandles.spectrum, 'YLim', [-140 -20]);
       set(PlotHandles.spectrum, 'XTick', [1, 2, 5, 10, 20, 50, 100]);
        set(PlotHandles.spectrum, 'XTickLabel', {'1.0', '2.0', '5.0', '10.0', '20.0', '50.0', '100.0'});
        xlabel('kHz');
        set(gca, 'XScale', 'log')
        
        rmsval = rms(ch2);
        MIC.RefSig
        MIC.Vref
        MIC
        rmsval
        dbspl = MIC.RefSig + 20.*log10(rmsval/MIC.Vref);
        % dbspl_spec = MIC.RefSig + 20*log10(2*Hpsd.Data/MIC.Vref); %#ok<NODEF>
        hdb = findobj('Tag', 'dbspl_text');
        length(dbspl)
        if(~isempty(hdb) && ishandle(hdb))
             set(hdb, 'String', sprintf('%7.1f', dbspl));
        end;

        tb = [0:1./AIsamp:1./AIsamp*(length(ch2)-1)];
        plot(PlotHandles.waveform, tb, ch2);
        
    case 'start' % --- Executes on button press in Start.
        dscale = 3600;
        if isempty(continueFlag)
            continueFlag = 0;
        end;
        if continueFlag == 0
            NoiseExposeNI('noisecalc');
            stdacrepeat = floor(Duration*3600.0); % duration in hours, convert to seconds
            elapsed = 0;
            totalDuration = Duration * dscale;
        else
            stdacrepeat = floor(remainingTime);
        end;
        armAO(stdacrepeat, swave, 0.);
        fprintf(2, 'total duration: %f   remaining time: %f\n', totalDuration, remainingTime);
        AO.Rate=AO_freq;
        if(size(swave, 1) > 1)
            AO.queueOutputData(swave);
        else
            AO.queueOutputData(swave');
        end;

    
        fprintf(1, 'Freq: %7.3f  BW: %8.3f  Attn: %3d\n', Freq, Band_Width, Attn);
        set_attn(Attn);
        armAO(stdacrepeat, swave, 0.);
%         AO.startForeground()
%         AO.TriggersPerRun=25;
%         addCounterInputChannel(AO,'Dev1','ctr0','EdgeCount');
        
%         a=get(AO);
%         if(isempty(a.Channel)) % only add if a channel does not already exist.
%             addchannel(AO, 1); % use dac 1 for output
%         end;
%         set(AO, 'samplerate', AO_freq); % STIM.NIFreq should is in samples per second.
%         set(AO, 'TriggerType', 'Immediate');
%         set(AO, 'RepeatOutput', stdacrepeat);
%         %  set(AO, 'BufferingMode', 'Auto');
%         %    set(AO, 'Timeout', 5*period);
%         if(size(swave, 1) > 1)
%             putdata(AO, swave);
%         else
%             putdata(AO, swave');
%         end;

        hdurshow = findobj('Tag', 'NoiseExposeTime');
        if(~isempty(hdurshow))
            set(hdurshow, 'String', sprintf('%02d:%02d:%02d   %6.1f   %5.1f%%', 0., 0., 0., 0., 0.));
        end;
        
        remainingTime = totalDuration;
        hdurleftshow = findobj('Tag', 'NoiseExposeRemaining');
        [rhr, rmins, rsec] = timecvt(remainingTime);
        if(~isempty(hdurleftshow))
            set(hdurleftshow, 'String', sprintf('%02d:%02d:%02d   %6.3f   %5.1f%%', rhr, rmins, rsec, ...
                remainingTime/3600., 100.*remainingTime/totalDuration));
        end;
        
%         AO.startForeground()
        
        lh = addlistener(AO,'DataRequired',@(src,event) src.queueOutputData(swave'));
        AO.IsContinuous = true;
        queueOutputData(AO,repmat(swave',5,1));
        AO.startBackground()
       % AO.wait()
        %         addCounterInputChannel(AO,'Dev1','ctr0','EdgeCount');

        %         start(AO); % get NI board read to go, then trigger the RP
        tstart = tic;
        QUIT = 0;
        continueFlag = 0; % reset the flag
        
        while QUIT == 0
            
            pause (1);
            elapsed = toc(tstart)+total_elapsed; % elapsed is in seconds into this overall epoch
            [hr, mins, sec ] = timecvt(elapsed);
            if(~isempty(hdurshow))
                set(hdurshow, 'String', sprintf('%02d:%02d:%02d   %6.3f   %5.1f%%', hr, mins, sec,...
                    elapsed/3600., 100.*elapsed/totalDuration));
            end;
            remainingTime = totalDuration-elapsed;
            [rhr, rmins, rsec] = timecvt(remainingTime);
            if(~isempty(hdurleftshow))
                set(hdurleftshow, 'String', sprintf('%02d:%02d:%02d   %6.3f   %5.1f%%', rhr, rmins, rsec, ...
                    remainingTime/3600., 100.*remainingTime/totalDuration));
            end;
            if(elapsed >= totalDuration) % convert to seconds
                QUIT = 1;
            end
        end
        set_attn(120);
        AO.stop()
        delete(lh)
        AO.release()
        %      invoke(RP, 'Halt');
        
    case 'stop' %--------------------------------
        
        set_attn(120.0); % set attenduation down
        stop(AO);
        %        invoke(RP, 'Halt');
        QUIT = 1;
        fprintf(1, 'Stopping!\n');
        return;
        
    case 'continue' % If there is time left on the elapsed clock from the last calculation
        % we pick up from where we ended... same as "start",
        % except the time comes from the remaining time...
        total_elapsed = elapsed;
        if elapsed < totalDuration
            continueFlag = 1;
            NoiseExposeNI('start');
        end;
    otherwise
end;
end

%-------------------------------------------------------------------
% Hardware interaction routines

% set_attn controls the PA5 Programmable attenuator
%
function set_attn(attn)
global PA5
global UseAttn

if(attn > 120 || attn < 0) % out of range : set to maximum attenuation
    attn = 120;
end;
if strncmp('left', UseAttn, 4)
    attn_id = 1;
else
    attn_id = 2;
end;
if(isempty(PA5))
    PA5=actxcontrol('PA5.x', [5 5 26 26]);
    if(PA5.ConnectPA5('USB', attn_id) == 0) % 1 is left, 2 is right
        fprintf(2, ' ... failed to connect to PA5\n');
        return;
    end;
    %     if(invoke(PA5, 'ConnectPA5', 'USB', 2) == 0)
    %         disp('failed to connect to PA5');
    %     end;
end;
invoke(PA5, 'SetAtten', attn);

hattn = findobj('Tag', 'Show_Atten');
if(~isempty(hattn))
    set(hattn, 'String', sprintf('%3d', attn));
end;

end

function [hr, mins, sec] = timecvt(elapsed)
hr = floor(elapsed/3600); % converts to hours
mins = mod(floor(elapsed-mod(elapsed,60))/60, 60); % minutes into current hour
sec = floor(mod(elapsed, 60));

end

function armAO(dursecs, swave, fs)
global AO Freq Band_Width Attn AO_freq

%set(AO, 'RepeatOutput', dursecs); TFR
%  set(AO, 'BufferingMode', 'Auto');
%    set(AO, 'Timeout', 5*period);
if(size(swave, 1) > 1)
    AO.queueOutputData(swave);
    %putdata(AO, swave);
else
    AO.queueOutputData(swave');
    %putdata(AO, swave');
end;

    
fprintf(1, 'Freq: %7.3f  BW: %8.3f  Attn: %3d\n', Freq, Band_Width, Attn);
set_attn(Attn);
end


function [ch1, ch2, AISamp, err] = checklevel(swave, fs)
%% Load and record signal to do a calibration
%% assumes stimulus waveform is already in swave.
%% sets up AO and RP .
%% 5/1/2010 Paul B. Manis
%%
global AO RP MicIn
% use a 2 second burst.
stimdur = 2.0; % stimulus duration in seconds
stop(AO);
armAO(stimdur, swave, fs)
nRecordPoints = 100000.;

if isempty(RP)
    RP=actxcontrol('rpco.x', [5 5 26 26]);
    if(invoke(RP, 'connectrp2', 'usb', 1) == 0)
        error('Failed to connect to RP2.1');
    end;
end

err = rp_setup(nRecordPoints);
AISamp = MicIn.sample_freq;
startBackground(AO); % get ni board read to go, then trigger the rp
pause(0.1);
RP.SoftTrg(1); % start.
curindex=RP.GetTagVal('Index');
while(curindex  < nRecordPoints) % Checks to see if it has read into half the buffer
    curindex=RP.GetTagVal('Index');
end
fprintf(1, 'RP done\n');
RP.SoftTrg(2);

ch1=double(RP.ReadTagV('data_out1', 0, nRecordPoints));
ch2=double(RP.ReadTagV('data_out2', 0, nRecordPoints));
stop(AO);
set_attn(-1);
RP.SoftTrg(2);
RP.Halt;
end

function setup_NI()
global AO 
daq.getDevices();
AO = daq.createSession('ni');
addAnalogOutputChannel(AO, 'Dev1', 1, 'Voltage');

% addDigitalChannel(NI, 'Dev1', 'Port0/Line5', 'InputOnly');
% PARS.ExposureChannel = 4;
% Counter = addCounterOutputChannel(AO, 'Dev1', 0, 'PulseGeneration');
AO.Rate = 500000; % Hz
Counter.Frequency = AO.Rate;

end

function err = rp_setup(nRecordPoints, varargin)
% TDT RP2.1 configuration routine
%
% 5/3/2010 Paul B. Manis
% Assumes that the RP2.1 had already been instantiated
% activeX

global RP STIM STOP MicIn
RP=actxcontrol('rpco.x', [5 5 26 26]);
if(invoke(RP, 'connectrp2', 'usb', 1) == 0)
    error('failed to connect to rp2');
end;
STIM.NIFreq = 500000;

% designate the sample frequency
% 0 = 6K, 1 = 12K, 2 = 25k, 3 = 50k, 4 = 100k, 5 = 200k, > 5 is not defined.
samp_cof_flag = 5; % 4 is for 100 kHz
samp_flist = [6103.5256125, 122107.03125, 24414.0625, 48828.125, ...
    97656.25, 195312.5];
err = 0;
MicIn.RP2COFFlag = samp_cof_flag;
MicIn.sample_freq = samp_flist(samp_cof_flag);
%fprintf(1, 'true sample frequency: %.6f hz\n', sfreq);

% see if the desired action is to STOP the RP2.1 from running
if nargin > 1 && strcmp(varargin{1}, 'stop')
    check_stop(1); % general STOP routine
    err = 1;
    return;
end

if RP.ClearCOF() == 0
    error('failed to clear cof');
end;
thisdir = pwd;
if (RP.LoadCOFsf([thisdir '\abrs.rcx'], MicIn.RP2COFFlag) == 0)
    error ('failed to load abrs.rcx file');
end;
RP.Run();
status = double(RP.GetStatus());
if bitget(double(status), 1) == 0;
    fprintf(2, 'rp_setup: Error connecting to RP2.1\n');
    err = 1;
    return;
elseif bitget(double(status), 2) == 0;
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
if RP.SetTagVal('REC_Size', nRecordPoints*2) == 0
    fprintf(2, 'rp_setup: Failed to set tag value rec_size with n_record_points\n');
    err = 1;
    return;
end;

if nargin > 1 && strcmp(varargin{1}, 'start')
    if RP.SoftTrg(1) == 0 % start.
        fprintf(2, 'failed to set trigger on RP2.1\n');
        set(hstat, 'string', 'rp error');
        check_stop(1);
        STOP = 0;
        err = 1;
        return;
    end;
end

end

