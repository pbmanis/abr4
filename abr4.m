function []=abr4(varargin)
% Program for Auditory Brainstem Evoked Response Measurements
% Paul B. Manis, Ph.D. UNC CHapel Hill, Otolaryngology/Head and Neck
% Surgery
% Updated: 2 March 2007
%
% Main code to collect abr (Auditory Brainstem Response) data
% Utilizing Tucker-Davis system 3 equipment and Active-X interface
% via USB port.
% Requires 14-bit NI board (PCI 6731) for output signal generation.

% OLD Calibration information:
% Etymotic ER7C, SN 87491, fresh battery.
% Output 94 dB SPL at 1 kHz. Measured microphone output directly.
% Reference tone yields 160 mV P-P, or 59.0 mV RMS
% dB = 20 * log(V/Vref)
% or:
% Vref = V * 10*(-dB/20)
% Vref = 0.0011772 V (0 dB SPL reference)
% using RMS voltage measurements.
%
% The click generates 540 mV P-P at 20 dB attenuation
% Referencing this to RMS


global STOP IN_ACQ STIM
global DATAa DATAp DATAn REFERENCE
global DLINE SIGNAL
global SCALE
global RESPONSE_MAP PLOTHANDLES
global CDATAp CDATAn
global ACQ4_FIG
global SPKR
global AO
global RP
global CAL SPLCAL %#ok<NUSED>
global HARDWARE
global DataDirectory


if (exist('last_data', 'var') == 0)
    last_data = 'None';
end
if (exist('GrandDatap', 'var') == 0)
    GrandDatap = {};
    GrandDatan = {};
end
DataDirectory = 'C:\Users\experimenters\Desktop\ABR_Data';

if(nargin == 0) % initialize
    clear HARDWARE
    clear STIM
    
    HARDWARE.system = computer(); % get the computer - this parses a few things.
    STIM.Info = 'ABR4 StimFile';
    STIM.StimPerSweep = 1;
    PA5 = [];
    ACQ4_FIG=open('ABR4.fig'); % get the figure window
    datacursormode();

    hardware_initialization(); % init the hardware.
%    DataDirectory = 'C:\Users\Experimenters\Desktop\ABR_Data';
    IN_ACQ = 0; % flag for when we are busy in acquisition
    DLINE = [];
    hf = [];
    RESPONSE_MAP = [];
    SIGNAL = [];
    SCALE = 0;
    [Speaker, Mic] = getSpeakerMic();
    SPKR.id = Speaker;
    SPKR.attn = 0;
    hm = findobj('Tag', 'ABR_cursor');
    msg_pos = get(hm, 'Position');
    set(gcf, 'UserData', []);
    PLOTHANDLES.signal1 = findobj('tag', 'ABR_Stimulus1');
    PLOTHANDLES.signal2 = findobj('tag', 'ABR_Stimulus2');
    PLOTHANDLES.data = findobj('tag', 'ABR_AvgData');
    PLOTHANDLES.responsemap = findobj('tag', 'ABR_ResponseMap');
    
    %     datac_setcrosshair(hf(1), 'ABR_stim1', 'ms', 'V', msg_pos);
    %     hf(2) = get_axis('signal2');
    %         datac_setcrosshair(hf(2), 'ABR_stim2', 'ms', 'V', msg_pos);
    %     hf(3) = get_axis('data', 1);
    %
    %     datac_setcrosshair(hf(3), 'ABR_data', 'ms', '\microV', msg_pos);
    %     hf(4) = get_axis('response');
    %
    %     datac_setcrosshair(hf(4), 'ABR_map', 'F (kHz)', 'SPL(dB)', msg_pos);
    %     set(gcf, 'UserData', hf); % store handles... use get_axis routine to recover handles
    set_attn(-1);
    dac_zero;
    STOP = 0;
    last_data = 'clickabr';
    return;
end

% always check the speaker prior to stimulation
[Speaker, ~] = getSpeakerMic();
SPKR.id = Speaker;
switch (SPKR.id(1:3))
    case {'ES1', 'EC1'}
        SPKR.attn = 0.0;
        SPLCAL.maxtones = 83.9; % New calibration, 5/1/2010 P. Manis Assumes ES Driver at - 6dB for linearity
        SPLCAL.maxclick = 79.5; % 84.8; 79 is with 6db attenuation ES1...
    case {'MF1', 'FF1'}
        SPKR.attn = 30.0; % for tones... 
        SPLCAL.maxtones = 110.0; % for mf1 speaker
        SPLCAL.maxclick = 108.5; % set with peak 1/4" mic output to match 80dB spl tone at "1e-6"
        % 114.8; % Old calibration 2007-4/30/2010db SPL with 0 dB attenuation (5 V signal)
    otherwise
        fprintf(2, 'Speaker type not known: %s %s\n', SPKR.id(1:3), SPKR.id);
        return;
end

% Handle the gui callbacks.
cmd = varargin{1};
switch(cmd)
    case 'quit' % do a clean shutdown
        close(ACQ4_FIG);
        clear;
        return;
    case 'load' % load a stimulus file from disk
        [FileName,PathName,~] = uigetfile('*.abr4','Stim File to Load', 'StimFiles/*.abr4');
        if FileName == 0 % cancelled out.
            return;
        end
        
        s = load([PathName FileName], '-mat');
        isfield(s.STIM, 'Info');
        s.STIM.Info
        if isfield(s.STIM, 'Info') && strcmp(s.STIM.Info, 'ABR4 StimFile')
            STIM = s.STIM;
        else
            fprintf(2, 'File does not appear to be an ABR4 Stim File\n');
            return;
        end
        updateStimParams(); % store info back to the window...
        hfn = findobj('tag', 'ABR_StimFile');
        if ~isempty(hfn)
            set(hfn, 'String', FileName);
        end
        
        
    case 'save' % save a stimulus file to disk
        [FileName,PathName] = uiputfile('*.abr4','Stim File to Save',...
            'StimFiles/*.abr4');
        if FileName == 0 % cancelled out.
            return;
        end
        getStimParams(); % read the current data in the window.
        save([PathName FileName], 'STIM');
        hfn = findobj('tag', 'ABR_StimFile');
        if ~isempty(hfn)
            set(hfn, 'String', FileName);
        end
        
    case {'calibrate', 'microphone', 'microphone104', 'checkcal'}
        % access calibration routines.
        acquire4(cmd);
        
    case 'response_spec' % calculate the response spectrum
        if(~isempty(DATAp))
            Hs=spectrum.periodogram('blackman');
            d = DATAp(~isnan(DATAp));
            Hpsd=psd(Hs,d, 'Fs', STIM.sample_freq/1000);
            % clear the axes of the response map-- TFR 11/16/2015
            cla(PLOTHANDLES.responsemap);
            plot(PLOTHANDLES.responsemap, Hpsd.Frequencies(2:end), ...
                log10(Hpsd.Data(2:end)));
            set(PLOTHANDLES.responsemap, 'XLim', [0.2, 64.]);
            set(PLOTHANDLES.responsemap, 'XScale', 'log');
            set(PLOTHANDLES.responsemap, 'YScale', 'linear');
            set(PLOTHANDLES.responsemap, 'YLimMode', 'auto');
            xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
            set(PLOTHANDLES.responsemap, 'XTick', xt);
            set(PLOTHANDLES.responsemap, 'XTickLabel', xt);
            
        end
        
    case 'stim_spec' % stimulus spectrum, ref max response
        if(~isempty(STIM.wave))
            Hs=spectrum.periodogram('blackman');
            Hpsd=psd(Hs,STIM.wave,'Fs', STIM.NIFreq/10);
            cla(PLOTHANDLES.responsemap);
            normSpec = Hpsd.Data(2:end)/max(Hpsd.Data(2:end));
            plot(PLOTHANDLES.responsemap, Hpsd.Frequencies(2:end), ...
                log10(normSpec));
            set(PLOTHANDLES.responsemap, 'XLim', [0.2, 64.0]);
            set(PLOTHANDLES.responsemap, 'XScale', 'log');
            set(PLOTHANDLES.responsemap, 'YScale', 'linear');
            set(PLOTHANDLES.responsemap, 'YLimMode', 'auto');
            xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
            set(PLOTHANDLES.responsemap, 'XTick', xt);
            set(PLOTHANDLES.responsemap, 'XTickLabel', xt);
        end
        
    case 'map' % redraw the response map
        cla(PLOTHANDLES.responsemap);
        switch (last_data)
            case 'clickabr'
                s = size(CDATAp);
                maxr = zeros(s(1),1);
                at=zeros(s(1), 1);
                for i = 1:s(1)
                    maxr(i) = max(abs(CDATAp(i,:)'+CDATAn(i,:)'));
                    at(i) = STIM.spls(i);
                end
                get_axis('response');
                plot(PLOTHANDLES.responsemap, at,  maxr,...
                    'marker', 'x', 'color', 'green');
                
            case 'toneabr'
                if(isempty(STIM.spls) || isempty(STIM.freqs))
                    return;
                end
                clear_plots;
                fr = STIM.freqs{1}(:);
                maxr = zeros(length(STIM.spls), length(fr));
                set(PLOTHANDLES.responsemap, 'XLim', [min(0.8*fr), max(1.2*fr)]);
                set(PLOTHANDLES.responsemap, 'YLim', [min(STIM.spls), max(STIM.spls)]);
                set(gca, 'Xscale', 'log');
                RESPONSE_MAP = quiver(fr, STIM.spls, 0*maxr, maxr); % set(RESPONSE_MAP, 'Xdata', fr, 'Ydata', spls, 'Zdata', maxr);
                set(RESPONSE_MAP, 'marker', '^', 'markersize', 1.5); % set(RESPONSE_MAP, 'Zdata', maxr);
                drawnow;
                for j = 1:length(fr) % go over the response map, frequency x atten
                    CDATAp = GrandDatap{j};
                    CDATAn = GrandDatan{j};
                    for i = 1:length(STIM.spls)
                        maxr(i,j) = [maxr max(abs(CDATAp(i,:)'+CDATAn(i,:)'))]; % measure the signal
                        get_axis('response');
                        r = quiver(fr, -STIM.spls, 0*maxr, maxr, 0.67);
                        set(r, 'marker', 'x', 'markersize', 1.5); % set(RESPONSE_MAP, 'Zdata', maxr);
                        drawnow
                    end
                end
                
                
            otherwise
        end
        
    case {'click', 'click_test'}  % acquire click abr series - intensity series
        if(IN_ACQ)
            return;
        end
        getStimParams;
        clear_plots;
        REFERENCE = [];
        hstat = findobj('tag', 'ABR_Status');
        cla(PLOTHANDLES.responsemap);
        set(PLOTHANDLES.responsemap, 'XScale', 'linear');
        set(PLOTHANDLES.responsemap, 'YLimMode', 'manual');
        set(PLOTHANDLES.responsemap, 'XLimMode', 'manual');
        set(PLOTHANDLES.responsemap, 'XLim', [min(STIM.spls), max(STIM.spls)]);
        set(PLOTHANDLES.responsemap, 'YLim', [0, 1]);
        STOP = 0;
        c = clock;
        fnamep = sprintf('%s/%4d%02d%02d-%02d%02d-p.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5));
        fnamen = sprintf('%s/%4d%02d%02d-%02d%02d-n.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5));
        
        fnamei = sprintf('%s/%4d%02d%02d-%02d%02d-SPL.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5));
        hf = findobj('Tag', 'ABR_CurrFreq');
        set(hf, 'String', 'Click');
        
        if strcmp(cmd, 'click_test')
            nspl = 1;
            spllist = [75];
            mode = 'test';
            hfn = findobj('Tag', 'ABR_filename');
            if(ishandle(hfn))
                set(hfn, 'string', '   ');
            end
        else
            nspl = length(STIM.spls);
            spllist = STIM.spls;
            mode = 'real';
            hfn = findobj('Tag', 'ABR_filename');
            if(ishandle(hfn))
                set(hfn, 'string', fnamep);
            end
        end
        at = NaN*ones(nspl, 1);
        maxr = NaN*zeros(nspl, 1);
        for i = 1 : nspl
            hf = findobj('Tag', 'ABR_CurrSPL');
            set(hf, 'String', sprintf('%.1f dB', spllist(end-(i-1))));
            stf = check_stop(0);
            if stf == 1 % successful STOP from the button
                err = 1;
                return;
            end
            err = click_abr(spllist(end-(i-1)), mode);
            if i == 1
                CDATAp = zeros(nspl, length(DATAp)); % positive and negative data arrays for click data
                CDATAn = zeros(nspl, length(DATAn));
            end    
            if(err > 0)
                return;
            end
            CDATAp(i,:) = DATAp; % append the averaged positive data
            CDATAn(i,:) = DATAn; % and the averaged negative data
            maxr(i) = max(abs(CDATAp(i,:)'+CDATAn(i,:)')); % measure the signal
            at(i) = spllist(end-(i-1));

            cla(PLOTHANDLES.responsemap);
            %tessa editing
            if ~isempty(maxr(~isnan(maxr)))
                plot(PLOTHANDLES.responsemap, at(~isnan(at)), ...
                maxr(~isnan(maxr))*1e6, 'ko-', ...
                'linewidth', 2, ...
                'markerfacecolor', 'red');
                set(PLOTHANDLES.responsemap, 'YLim', [0, max(maxr(~isnan(maxr))*1e6)]);
            else 
                plot(PLOTHANDLES.responsemap, at(~isnan(at)), ...
                0, 'ko-', ...
                'linewidth', 2, ...
                'markerfacecolor', 'red');
               
            end
            
            
        end
        set(hstat, 'String', 'Done');
        CDATAp = CDATAp';
        CDATAn = CDATAn';
        spl = STIM.spls';
        if strcmp(cmd, 'click')
            save(fnamep, 'CDATAp', '-ascii', '-tabs');
            save(fnamei, 'spl', '-ascii', '-tabs'); % save intensity list also
            save(fnamen, 'CDATAn', '-ascii', '-tabs');
        end
        last_data = cmd;
        
        
        % all TONE stimuli(F, I maps or single points) go through here...
    case {'tone_abr', 'tone_mapping', 'tone_test', 'tone_info'} % execute the tone ABR measurement
        if(IN_ACQ)
            return;
        end
        clear_plots;
        STOP = 0;
        getStimParams;
        % Get the calibration frequency map for the current speaker
        %
        load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file. Result is in CAL
        if(isempty(STIM.spls) || isempty(STIM.freqs))
            return;
        end
        hstat = findobj('tag', 'ABR_Status');
        
        if strcmp(cmd, 'tone_test')
            spllist = [75]; %#ok<*NBRAK>
            fr = [8000];
            mode = 'test';
        elseif strcmp(cmd, 'tone_info')
            nspl = length(STIM.spls);
            spllist = STIM.spls;
            fr = STIM.freqs{1}(:);
            mode = 'info';
        else
            nspl = length(STIM.spls);
            spllist = STIM.spls;
            fr = STIM.freqs{1}(:);
            mode = 'real';
        end
        GrandDatap = cell(length(fr), 1);
        GrandDatan = cell(length(fr), 1);
        
        maxr = zeros(length(STIM.spls), length(fr));
        maxry = maxr;
        maxrx = maxr;
        if ~strcmp(cmd, 'tone_test')
            hr = PLOTHANDLES.responsemap;
            cla(hr);
            set(hr, 'YLimMode', 'auto');
            set(hr, 'XLimMode', 'manual');
            %            set(hr, 'XLim', [min(0.8*fr), max(1.2*fr)]);
            set(hr, 'YLim', [min(STIM.spls)-5, max(STIM.spls)]+5);
            set(PLOTHANDLES.responsemap, 'XLim', [0.2, 64.]);
            set(PLOTHANDLES.responsemap, 'XScale', 'log');
            xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
            set(PLOTHANDLES.responsemap, 'XTick', xt);
            set(PLOTHANDLES.responsemap, 'XTickLabel', xt);
            
            RESPONSE_MAP = quiver(hr, fr/1000.0, STIM.spls, 0*maxr, maxr, 1e-3);
            set(RESPONSE_MAP, 'marker', 'o', 'markersize', 0.8,'markerfacecolor', 'k', 'markeredgecolor', 'k');
            set(RESPONSE_MAP, 'AutoScale', 'off');
            set(RESPONSE_MAP, 'UDataSource', 'maxrx');
            set(RESPONSE_MAP, 'VDataSource', 'maxry');
            refreshdata(RESPONSE_MAP, 'caller');
        end
        err = 0;
        c=clock;
        %  fprintf(1, 'nattn: %d,  nfreq: %d\n', length(spls), length(fr));
        fnamei = sprintf('%s/%4d%02d%02d-%02d%02d-SPL.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5));
        fnamef = sprintf('%s/%4d%02d%02d-%02d%02d-kHz.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5));
        hf = findobj('Tag', 'ABR_CurrFreq');
        hs = findobj('Tag', 'ABR_CurrSPL');
        hfn = findobj('Tag', 'ABR_filename');
        for j = 1:length(fr) % go over the response map, frequency x atten
            CDATAp = [];
            CDATAn = [];
            fnamep = sprintf('%s/%4d%02d%02d-%02d%02d-p-%8.3f.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5), fr(j));
            fnamen = sprintf('%s/%4d%02d%02d-%02d%02d-n-%8.3f.txt', DataDirectory, c(1), c(2), c(3), c(4), c(5), fr(j));
            set(hf, 'String', sprintf('%6.1f kHz', fr(j)));
            if(~strcmp( cmd, 'tone_test'))
                if(ishandle(hfn))
                    set(hfn, 'string', fnamep);
                end
            else
                set(hfn, 'String', 'Tone Test (no file)');
            end
            
            for i = 1:length(spllist)
                stf = check_stop(0);
                if stf == 1 % successful STOP from the button
                    err = 1;
                    return;
                end

                err = tone_map(cmd, fr(j), spllist(i)); % main tone abr routine
                if(err ~= 0)
                    return;
                end
                set(hs, 'String', sprintf('%3.1f dB', spllist(i)));
                
                if ~strcmp(cmd, 'tone_test')
                    maxr(i,j) = max(abs(DATAa(1,:)')); % measure the signal
                    tmax = max(max(maxr));
                    maxry = 5*maxr/tmax;
                    %                fprintf(2, 'maxr(i,j): %g   tmax: %g\n', maxr(i,j), tmax);
                    refreshdata(RESPONSE_MAP, 'caller');
                    drawnow
                end
                CDATAp = [CDATAp; DATAp]; % append the positive data
                CDATAn = [CDATAn; DATAn]; % and the negative data
                %   fprintf(1, 'fr = %6.1f, j = %d, spl = %6.1f  i = %d\n', fr(j), j, spls(i), i);
                %               pause(0.1); % brief delay between frequencies.
            end
            GrandDatap{j} = CDATAp;
            GrandDatan{j} = CDATAn;
            CDATAp = CDATAp';
            CDATAn = CDATAn';
            if(~strcmp(cmd, 'tone_test'))
                save(fnamep, 'CDATAp', '-ascii', '-tabs');
                save(fnamen, 'CDATAn', '-ascii', '-tabs');
            end
            
        end
        spl = spllist';
        if ~strcmp(cmd, 'tone_test')
            save(fnamei, 'spl', '-ascii', '-tabs'); % save intensity list also
            save(fnamef, 'fr', '-ascii', '-tabs'); % and frequency list
        end
        set(hstat, 'String', 'Done');
        last_data = cmd;
        
        
    case 'abort'
        STOP = 1;
        last_data = 'None'; %#ok<*NASGU>
        fprintf(2, 'Abort hit\n');
    otherwise
end
end


%----------------------------------------------------------------------------
%****************************************************************************
%----------------------------------------------------------------------------

% click_abr
% function to generate clicks and collect abr results
% This routine collects data for one stimulus condition (SPL)
% 9/23/03 P. Manis
%

function [err] = click_abr(spl, mode)
% take a single parameter (spl, mode) click abr
% requires call to getStimParams first.
% if mode is 'test', we run in a reduced acquisition format
%
global STIM DATAa DATAp DATAn REFERENCE SPLCAL AO HARDWARE PLOTHANDLES

fprintf(2, 'Calling click_abr\n');
err = 0;
STIM.rate = 1/STIM.sample_freq; % rate is in sec per point (recording side).
STIM.click_dur = 0.1;  % default is 50 microseconds (0.05)
keep_reference = 1;
if ~exist('REFERENCE', 'var')
    REFERENCE = [];
end

if strcmp(mode, 'test')
    oldsweeps = STIM.NSweeps;
    STIM.NSweeps = 10;
    oldsps = STIM.StimPerSweep;
    STIM.StimPerSweep = 40;
    oldipi = STIM.ipi;
    STIM.ipi = 25;
    set(STIM.hsr, 'String', num2str(STIM.NSweeps));
    updateStimParams;
end

cnp = STIM.StimPerSweep;

if strcmp(HARDWARE, 'NI')
    STIM.NIFreq = 500000; % express in sec per point
else
    STIM.NIFreq = 44100;
end
STIM.delay = STIM.click_delay;
updateStimParams;

[STIM.wave, STIM.clock] = click(SPLCAL.click_amp, STIM.click_delay, STIM.click_dur,...
    STIM.NIFreq, STIM.ipi, cnp, STIM.Alternate); % convert rate to usec per point
attn = SPLCAL.maxclick - spl;
if(attn < 0)
    fprintf(1, 'Requesting sound louder than available with this speaker\nSetting to 0 attn\n');
    attn = 0;
end
clear_plots(keep_reference);
fl = get_SignalPlotFlag(1);
STIM.Monitor = get_SignalPlotFlag(2);
if fl == 1
    ts = (0:length(STIM.wave)-1)*1000*STIM.clock;
    maxt = (length(STIM.wave)-1)*1000*STIM.clock; % in msec
%     figure; plot(ts,PLOTHANDLES.signal1);
    plot(PLOTHANDLES.signal1, ts, STIM.wave, 'color', 'blue');
    set(PLOTHANDLES.signal1, 'XLim', [0 maxt]);
    drawnow;
end

[data, err] = acquire4('attn', attn);
set_attn(-1);
if(err == 0)
    davg = final_plot(data, keep_reference);
    DATAp = data(1,:);
    DATAn = data(2,:);
    DATAa = davg;
else
    fprintf(2, 'acquire4 err: %d\n', err);
end
if keep_reference && err == 0
    REFERENCE = davg;
end
if strcmp(mode, 'test')
    STIM.NSweeps = oldsweeps;
    STIM.ipi = oldipi;
    STIM.StimPerSweep = oldsps;
    updateStimParams;
end

end


%----------------------------------------------------------------------------
%****************************************************************************
%----------------------------------------------------------------------------

function [err] = tone_map(mode, freq, spl)

global STIM DATAa DATAp DATAn REFERENCE SPLCAL AO
global  HARDWARE PLOTHANDLES CAL SPKR

%[Speaker, Mic] = getSpeakerMic();
err = 0;
STIM.rate = 1/STIM.sample_freq; % rate is in sec per point.
if strcmp(HARDWARE, 'NI')
    STIM.NIFreq = 500000; % express in sec per point
else
    STIM.NIFreq = 44100;
end

if ~isempty(AO)
    if(STIM.Alternate)
        tnp = STIM.StimPerSweep;
    else
        tnp = STIM.StimPerSweep;
    end
else
    tnp = STIM.StimPerSweep;
end
if strcmp(mode, 'tone_test')
    REFERENCE = [];
    oldsweeps = STIM.NSweeps;
    STIM.NSweeps = 5;
    oldsps = STIM.StimPerSweep;
    STIM.StimPerSweep = 40;
    oldipi = STIM.ipi;
    STIM.ipi = 25;
    set(STIM.hsr, 'String', num2str(STIM.NSweeps));
end

STIM.delay = STIM.tone_delay;
updateStimParams;
% interpolate to get the attenuation at the requested frequency
% splatF=interp1(CAL.Freqs, CAL.maxdB, freq, 'spline'); 
splatF = soundfuncs.spl_at_f(CAL.Freqs, CAL.maxdB, freq);
attn = splatF-spl+SPKR.attn;
if(attn <= 0)
    fprintf(1, 'Requesting sound louder than available with this speaker\nSetting to 0 attn\n');
    return;
end
fprintf(2, 'mode: %s', mode)
if strcmp(mode, 'tone_info')
    fprintf(2, 'Freq: %8.2f  SPL: %8.2f  Attn: %6.1f\n', freq, spl, attn);
    return;
end

switch (mode)
    case {'tone_abr', 'test', 'tone_test'}
        disp ' tone_abr'
        [STIM.wave, STIM.clock] = tonepip(SPLCAL.ToneVoltage, freq, STIM.tone_delay, ...
            5, 0.5, 0, STIM.NIFreq, STIM.ipi, tnp, STIM.Alternate); % convert rate to usec per point
    case {'tone_mapping', 'tone_info'}
        disp 'tone_mapping'
        htd = findobj('tag', 'ABR_ToneMapDur');
        sel = get(htd, 'value');
        slist = str2num(get(htd, 'String')); %#ok<ST2NM>
        tdur = slist(sel);
        [STIM.wave, STIM.clock] = tonepip(SPLCAL.ToneVoltage, freq, STIM.tone_delay_mapping, tdur, 2.5, 0, STIM.NIFreq, STIM.ipi, STIM.StimPerSweep, 0); % convert rate to usec per point
    otherwise
        return;
end


clear_plots;
fl = get_SignalPlotFlag(1);
STIM.Monitor = get_SignalPlotFlag(2);
if fl == 1
    ts = (0:length(STIM.wave)-1)*1000*STIM.clock;
    maxt = (length(STIM.wave)-1)*1000*STIM.clock; % in msec
    plot(PLOTHANDLES.signal1, ts, STIM.wave, 'color', 'blue');
    set(PLOTHANDLES.signal1, 'XLim', [0 maxt]);
    set(PLOTHANDLES.signal1, 'clipping', 'on');
    
    drawnow;
end

[data, err] = acquire4('attn', attn);
set_attn(-1);
if (err ~= 0)
    return
end
if(err == 0)
    davg = final_plot(data);
    DATAp = data(1,:);
    DATAn = data(2,:);
    DATAa = davg;
end
if(strcmp(mode, 'tone_test'))
    STIM.NSweeps = oldsweeps;
    STIM.ipi = oldipi;
    STIM.StimPerSweep = oldsps;
    updateStimParams;
end
end

function [davg] = final_plot(data, varargin)
%tr = (0:length(data)-1)*STIM.rate*1000; % express rate in usec
global STIM ACQPars PLOTHANDLES REFERENCE
if(STIM.Alternate)
    davg = (data(1,:)+data(2,:))/2.0; % /(ACQPars.np+ACQPars.nn);
    dsum = (data(1,:)-data(2,:))/2.0;
else
    davg = data(1,:); % already the average
    dsum = davg;
end
if nargin > 0
    REFERENCE = davg;
end
hold on;
plot(PLOTHANDLES.data,  ACQPars.tb, davg, 'k-');
%md = max(abs(dsum));
%trace_scale(md, PLOTHANDLES.data);
end







