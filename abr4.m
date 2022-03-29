function [D] = abr4(varargin)

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

    persistent DataDirectory
    persistent HW
    persistent CALIBRATION
    persistent PARS
    persistent PLOTS
    persistent DATA
    persistent STIM
    persistent GUI

    D = []; % in case ABR4 needs to return something.

    if isempty(varargin) || (exist('HW', 'var') == 0)
        HW = abr4_hardware_struct;
        HW.initialize();

        STIM = abr4_STIM_struct;
        STIM.initialize();

        CALIBRATION = abr4_calibration_struct;
        CALIBRATION.initialize();

        PARS = abr4_parameters_struct;
        PARS.initialize();

        PLOTS = abr4_plots_struct;
        PLOTS.initialize();

        DATA = abr4_data_struct;
        DATA.initialize();

        GUI = abr4_gui_struct; % don't initialize until the gui is built.

        % fprintf(2, 'Created structures\n');

    end

    if (exist('last_data', 'var') == 0)
        last_data = 'None';
    end

    if (exist('GrandDatap', 'var') == 0)
        GrandDatap = {};
        GrandDatan = {};
    end

    DataDirectory = 'C:\Users\experimenters\Desktop\ABR_Data';
    CalibrationDirectory = 'calibration_history';

    if (nargin == 0) % initialize when entered for the first time.

        if ~exist(CalibrationDirectory, 'dir')
            mkdir(CalibrationDirectory);
        end

        clear HW.HARDWARE

        HW.HARDWARE.system = computer(); % get the computer - this parses a few things.
        STIM.Info = 'ABR4 StimFile';
        STIM.StimPerSweep = 1;
        PA5 = [];
        PARS.ABR4_FIG = open('ABR4.fig'); % get the figure window
        datacursormode();

        [HW, STIM, CALIBRATION] = hardware_initialization(HW, STIM, CALIBRATION); % init the hardware.
        %    DataDirectory = 'C:\Users\Experimenters\Desktop\ABR_Data';
        HW.IN_ACQ = 0; % flag for when we are busy in acquisition

        hf = [];
        PLOTS.responsemap = [];
        [Speaker, Mic] = getSpeakerMic();
        CALIBRATION = get_calibration_info(Speaker, Mic, CALIBRATION);

        CALIBRATION.SPKR.id = Speaker;
        msg_pos = get(GUI.hcursor, 'Position');
        set(gcf, 'UserData', []);
        PLOTS.signal1 = findobj('tag', 'ABR_Stimulus1');
        PLOTS.signal2 = findobj('tag', 'ABR_Stimulus2');
        PLOTS.data = findobj('tag', 'ABR_AvgData');
        PLOTS.responsemap = findobj('tag', 'ABR_ResponseMap');
        %     datac_setcrosshair(PLOT.signal1, 'ABR_stim1', 'ms', 'V', msg_pos);
        %     datac_setcrosshair(PLOT.signal2, 'ABR_stim2', 'ms', 'V', msg_pos);
        %     datac_setcrosshair(PLOT.data, 'ABR_data', 'ms', '\microV', msg_pos);
        %     datac_setcrosshair(PLOT.responsemap, 'ABR_map', 'F (kHz)', 'SPL(dB)', msg_pos);
        HW = set_attn(HW, -1);
        dac_zero;
%         STOP = 0;
%         last_data = 'clickabr';
        GUI = GUI.initialize();
        STIM = getStimParams(STIM, GUI); % now we can populate the STIM from the GUI
        set(GUI.hstop, 'UserData', 'Waiting');
        return;
    end

    % always check the speaker immediately prior to stimulation
    [Speaker, ~] = getSpeakerMic();
    CALIBRATION.SPKR.id = Speaker;

    switch (CALIBRATION.SPKR.id)
        case {'ES1', 'EC1'}
            CALIBRATION.SPKR.CalAttn = 0.0;
            CALIBRATION.SPLCAL.maxtones = 83.9; % New calibration, 5/1/2010 P. Manis Assumes ES Driver at - 6dB for linearity
            CALIBRATION.SPLCAL.maxclick = 79.5; % 84.8; 79 is with 6db attenuation ES1...
        case {'MF1'}
            CALIBRATION.SPKR.CalAttn = 20.0; % for tones...
            CALIBRATION.SPLCAL.maxtones = 110.0; % for mf1 speaker
            CALIBRATION.SPLCAL.maxclick = 108.5; % set with peak 1/4" mic output to match 80dB spl tone at "1e-6"
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  26 January 2022 pbmanis Re-calibration of the ABR system. click
            %  calibrated by averaging the microphone waveform (click_cal.m)
            % with 10 dB attenuation. This yields a value of 90.0 dB SPL
            % (measured peak to peak), at 0 attenaution. Click p-p voltage 22.1
            % mV at 80.03 dB (10 dB atten). Tone (16 kHz) p-p voltage 10.37 mV
            % or 73.6 dB SPL (requested 75). Calibrations for tone done on the
            % same day (26 Jan).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CALIBRATION.SPLCAL.maxclick = 90.0; % 26 Jan 2022. Note 18.5 dB difference from prior.
            % 114.8; % Old calibration 2007-4/30/2010db SPL with 0 dB attenuation (5 V signal)
        otherwise
            fprintf(2, 'Speaker type not known\n');
            return;
    end

    % Handle the gui callbacks.
    cmd = varargin{1};
    fprintf(2, "cmd: %s\n", cmd);

    switch (cmd)
        case 'quit' % do a clean shutdown
            close(PARS.ABR4_FIG);
            clear;
            return;
            
        case 'gethw' % return the HARDWARE, STIM and CALIBRATION structures
            D = struct("HW", HW, "STIM", STIM, "CALIBRATION", CALIBRATION);
            return
            
        case 'Stop'
            state = get(GUI.hstop, 'UserData');
            if nargin >= 2
                fprintf('abr4: STOP state: %s, message: %s\n', state, varargin{2})
            else
                fprintf('abr4: STOP state: %s, no message\n', state)
            end
            if strcmp(state, 'Waiting')
                set_status('Waiting');
                return
            elseif strcmp(state, 'Stopped')
%                 HW = stop_acq(HW, 'Stop');
%                 set_status('Stopped');
                return
            elseif strcmp(state, 'Running')
                HW = stop_acq(HW, 'Stop');
                set_status('Stopped');
            elseif strcmp(state, 'Paused')
                HW = stop_acq(HW, 'Stop');
                set_status('Paused');
            end
            
        case 'load' % load a stimulus file from disk
            STIM = loadStim(STIM, GUI);

        case 'save' % save a stimulus file to disk
            [FileName, PathName] = uiputfile('*.abr4', 'Stim File to Save', ...
            'StimFiles/*.abr4');

            if FileName == 0 % cancelled out.
                return;
            end

            [STIM] = getStimParams(STIM, GUI); % read the current data in the window.
            save([PathName FileName], 'STIM');
            hfn = findobj('tag', 'ABR_StimFile');

            if ~isempty(GUI.hstimfilename)
                set(GUI.hstimfilename, 'String', FileName);
            end

        case {'calibrate', 'microphone', 'microphone104', 'checkcal'}
            % access calibration routines.
            acquire4(cmd, HW, STIM, PLOTS, GUI, CALIBRATION);

        case 'response_spec' % calculate the response spectrum

            if (~isempty(DATAp))
                %             Hs=spectrum.periodogram('blackman');
                [Hpsd, w] = periodogram(DATA.DATAp, 'blackman', 'Fs', STIM.NIFreq / 10);
                %             d = DATAp(~isnan(DATAp));
                %             Hpsd=psd(Hs,d, 'Fs',STIM.sample_freq/1000);
                % clear the axes of the response map-- TFR 11/16/2015
                cla(PLOTS.responsemap);
                plot(PLOTS.responsemap, w, log10(Hpsd));
                set(PLOTS.responsemap, 'XLim', [0.2, 64.]);
                set(PLOTS.responsemap, 'XScale', 'log');
                set(PLOTS.responsemap, 'YScale', 'linear');
                set(PLOTS.responsemap, 'YLimMode', 'auto');
                xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
                set(PLOTS.responsemap, 'XTick', xt);
                set(PLOTS.responsemap, 'XTickLabel', xt);

            end

        case 'stim_spec' % stimulus spectrum, ref max response

            if (~isempty(STIM.wave))
                %             Hs=spectrum.periodogram('blackman');
                [Hpsd, w] = periodogram(STIM.wave, 'blackman', 'Fs', STIM.NIFreq / 10);
                cla(PLOTS.responsemap);
                %             normSpec = Hpsd.Data(2:end)/max(Hpsd.Data(2:end));
                plot(PLOTS.responsemap, w, log10(Hpsd));
                set(PLOTS.responsemap, 'XLim', [0.2, 64.0]);
                set(PLOTS.responsemap, 'XScale', 'log');
                set(PLOTS.responsemap, 'YScale', 'linear');
                set(PLOTS.responsemap, 'YLimMode', 'auto');
                xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
                set(PLOTS.responsemap, 'XTick', xt);
                set(PLOTS.responsemap, 'XTickLabel', xt);
            end

        case 'map' % redraw the response map
            cla(PLOTS.responsemap);

            switch (last_data)
                case 'clickabr'
                    s = size(DATA.CDATAp);
                    maxr = zeros(s(1), 1);
                    at = zeros(s(1), 1);

                    for i = 1:s(1)
                        maxr(i) = max(abs(DATA.DATAp(i, :)' + DATA.CDATAn(i, :)'));
                        at(i) = STIM.spls(i);
                    end

                    % get_axis('response');
                    plot(PLOTS.responsemap, at, maxr, ...
                    'marker', 'x', 'color', 'green');

                case 'toneabr'

                    if (isempty(STIM.spls) || isempty(STIM.freqs))
                        return;
                    end

                    clear_plots(PLOTS, STIM, DATA);
                    fr = STIM.freqs{1}(:);
                    maxr = zeros(length(STIM.spls), length(fr));
                    set(PLOTS.responsemap, 'XLim', [min(0.8 * fr), max(1.2 * fr)]);
                    set(PLOTS.responsemap, 'YLim', [min(STIM.spls), max(STIM.spls)]);
                    set(gca, 'Xscale', 'log');
                    PLOTS.responsemap = quiver(fr, STIM.spls, 0 * maxr, maxr); % set(PLOTS.responsemap, 'Xdata', fr, 'Ydata', spls, 'Zdata', maxr);
                    set(PLOTS.responsemap, 'marker', '^', 'markersize', 1.5); % set(PLOTS.responsemap, 'Zdata', maxr);
                    drawnow;

                    for j = 1:length(fr) % go over the response map, frequency x atten
                        CDATAp = GrandDatap{j};
                        CDATAn = GrandDatan{j};

                        for i = 1:length(STIM.spls)
                            maxr(i, j) = [maxr max(abs(CDATAp(i, :)' + CDATAn(i, :)'))]; % measure the signal
                            % get_axis('response');
                            r = quiver(fr, -STIM.spls, 0 * maxr, maxr, 0.67);
                            set(r, 'marker', 'x', 'markersize', 1.5); % set(PLOTS.responsemap, 'Zdata', maxr);
                            drawnow
                        end

                    end

                otherwise
            end

        case {'click', 'click_test'} % acquire click abr series - intensity series

            if (HW.IN_ACQ)
                return;
            end

            [STIM] = getStimParams(STIM, GUI);
            clear_plots(PLOTS, STIM, DATA);
            REFERENCE = [];
            cla(PLOTS.responsemap);
%             set(PLOTS.responsemap, 'XScale', 'linear');
%             set(PLOTS.responsemap, 'YLimMode', 'manual');
%             set(PLOTS.responsemap, 'XLimMode', 'manual');
%            set(PLOTS.responsemap, 'XLim', [min(STIM.spls), max(STIM.spls)]);
%            set(PLOTS.responsemap, 'YLim', [0, 1]);

            c = clock;
            basefilename = sprintf('%4d%02d%02d-%02d%02d', c(1), c(2), c(3), c(4), c(5));
            set(GUI.hdatafilename, 'String', basefilename);
            fnamep = sprintf('%s/%s-p.txt', DataDirectory, basefilename);
            fnamen = sprintf('%s/%s-n.txt', DataDirectory, basefilename);
            fname_bigdata = sprintf('%s/%s.mat', DataDirectory, basefilename);
            fnamei = sprintf('%s/%s-SPL.txt', DataDirectory, basefilename);
            
            set(GUI.hcurrentfrequency, 'String', 'Click');

            if strcmp(cmd, 'click_test')
%                 spllist = [75, 75, 75, 75, 75, 75, 75, 75, 75, 75];
                spllist = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90];
                nspl = length(spllist);
                mode = 'test';

                if (ishandle(GUI.hstimfilename))
                    set(hstimfilename, 'string', '   ');
                end

            else
                nspl = length(STIM.spls);
                spllist = STIM.spls;
                mode = 'real';

                if (ishandle(GUI.hstimfilename))
                    set(GUI.hstimfilename, 'string', '   ');
                end

            end

            STIM.sound_type = 'click';
            attns = NaN * ones(nspl, 1);
            maxr = NaN * zeros(nspl, 1);

            for i = 1:nspl
                set(GUI.hcurrentSPL, 'String', sprintf('%.1f dB', spllist(end - (i - 1))));
                state = get(GUI.hstop, 'UserData');
                abr4(state, 'click gui stop');
                if strcmp(state, 'Stop')
                    HW = stop_acq(HW, state);
                    set_status('Stopped');
                    err = 1;
                    return;
                end

                [DATA, STIM, err] = one_click(spllist(end - (i - 1)), mode, HW, STIM, CALIBRATION, DATA, PLOTS, GUI);

                if i == 1
                    DATA.CDATAp = zeros(nspl, length(DATA.DATAp)); % positive and negative data arrays for click data
                    DATA.CDATAn = zeros(nspl, length(DATA.DATAn));
                    DATA.C_CHDATA = zeros(nspl, length(DATA.CHDATA));
                end

                if (err > 0)
                    return;
                end

                DATA.CDATAp(i, :) = DATA.DATAp; % append the averaged positive data
                DATA.CDATAn(i, :) = DATA.DATAn; % and the averaged negative data
                DATA.C_CHDATA(i, :) = DATA.CHDATA; % raw traces 
                maxr(i) = max(abs(DATA.CDATAp(i, :)' + DATA.CDATAn(i, :)')); % measure the signal
                attns(i) = spllist(end - (i - 1));

                cla(PLOTS.responsemap);
                %tessa editing
                if ~isempty(maxr(~isnan(maxr)))
                    plot(PLOTS.responsemap, attns(~isnan(attns)), ...
                        maxr(~isnan(maxr)) * 1e6, 'ko-', ...
                        'linewidth', 2, ...
                        'markerfacecolor', 'red');
                    set(PLOTS.responsemap, 'YLim', [0, max(maxr(~isnan(maxr)) * 1e6)]);
                else
                    plot(PLOTS.responsemap, attns(~isnan(attns)), ...
                        0, 'ko-', ...
                        'linewidth', 2, ...
                        'markerfacecolor', 'red');

                end
                state = get(GUI.hstop, 'UserData');
                if strcmp(state, 'Stop')
                    HW = stop_acq(HW, state);
                    set_status('Stopped');
                    err = 1;
                    return;
                end
            end  % of the SPL loop for clicks

            HW = stop_acq(HW, state);
            % set(GUI.hstat, 'String', 'Done');
            DATA.CDATAp = DATA.CDATAp';
            DATA.CDATAn = DATA.CDATAn';
          %  DATA.C_CHDATA = DATA.C_CHDATA';
            spl = STIM.spls';

            if strcmp(cmd, 'click')
                CDATAp = DATA.CDATAp; % can't save part of a structure
                CDATAn = DATA.CDATAn;
                save(fnamep, 'CDATAp', '-ascii', '-tabs');
                save(fnamei, 'spl', '-ascii', '-tabs'); % save intensity list also
                save(fnamen, 'CDATAn', '-ascii', '-tabs');
                % make a structure with all the data and parameters
                bigdata = struct('CAL', CALIBRATION, 'DATA', DATA, 'STIM', STIM, 'HARDWARE', HW);
                save(fname_bigdata, 'bigdata', '-mat');
            end

            last_data = cmd;

        % all TONE stimuli(F, I maps or single points) go through here...
        case {'tone_abr', 'tone_mapping', 'tone_test', 'tone_info'} % execute the tone ABR measurement
            
            if HW.IN_ACQ
                return;
            end
% 
%             if CALIBRATION.Needs_Cal
%                 return
%             end

            clear_plots(PLOTS, STIM, DATA);
            [STIM] = getStimParams(STIM, GUI);

            if (isempty(STIM.spls) || isempty(STIM.FreqList))
                return;
            end

            hstat = findobj('tag', 'ABR_Status');

            if strcmp(cmd, 'tone_test')
                spllist = [90, 90, 90, 90, 90, 90, 90, 90, 90, 90];
                fr = [8000, 8000, 8000, 8000, 8000, 8000, 8000, 8000, 8000, 8000];
                mode = 'test';
            elseif strcmp(cmd, 'tone_info')
                nspl = length(STIM.spls);
                spllist = STIM.spls;
                fr = STIM.freqs;
                mode = 'info';
            else
                nspl = length(STIM.spls);
                spllist = STIM.spls;
                fr = STIM.freqs;
                mode = 'real';
            end

            GrandDatap = cell(length(fr), 1);
            GrandDatan = cell(length(fr), 1);

            maxr = zeros(length(STIM.spls), length(fr));
            maxry = maxr;
            maxrx = maxr;

            if ~strcmp(cmd, 'tone_test')
                hr = PLOTS.responsemap;
                cla(hr);
                set(hr, 'YLimMode', 'auto');
                set(hr, 'XLimMode', 'manual');
                %            set(hr, 'XLim', [min(0.8*fr), max(1.2*fr)]);
%                set(hr, 'YLim', [min(STIM.spls) - 5, max(STIM.spls)] + 5);
                set(PLOTS.responsemap, 'XLim', [0.2, 64.]);
                set(PLOTS.responsemap, 'XScale', 'log');
                xt = [0.5, 1, 2, 4, 8, 16, 32, 64];
                set(PLOTS.responsemap, 'XTick', xt);
                set(PLOTS.responsemap, 'XTickLabel', xt);

                PLOTS.quiver = quiver(hr, fr / 1000.0, STIM.spls, 0 * maxr, maxr, 1e-3);
                set(PLOTS.quiver, 'marker', 'o', 'markersize', 0.8, 'markerfacecolor', 'k', 'markeredgecolor', 'k');
                set(PLOTS.quiver, 'AutoScale', 'off');
                set(PLOTS.quiver, 'UDataSource', 'maxrx');
                set(PLOTS.quiver, 'VDataSource', 'maxry');
               % refreshdata(PLOTS.responsemap, 'caller');
            end

            err = 0;
            c = clock;
            %  fprintf(1, 'nattn: %d,  nfreq: %d\n', length(spls), length(fr));
            c = clock;
            basefilename = sprintf('%4d%02d%02d-%02d%02d', c(1), c(2), c(3), c(4), c(5));
            set(GUI.hdatafilename, 'String', basefilename);

            fnamei = sprintf('%s/%s-SPL.txt', DataDirectory, basefilename);
            fnamef = sprintf('%s/%s-kHz.txt', DataDirectory, basefilename);

            for j = 1:length(fr) % go over the response map, frequency x atten
                DATA.CDATAp = [];
                DATA.CDATAn = [];
                DATA.C_CHDATA = [];
                fnamep = sprintf('%s/%s-p-%8.3f.txt', DataDirectory, basefilename, fr(j));
                fnamen = sprintf('%s/%s-n-%8.3f.txt', DataDirectory, basefilename, fr(j));
                fname_bigdata = sprintf('%s/%s-%.1f.mat', DataDirectory, basefilename, fr(j));
                set(GUI.hdatafilename, 'String', fname_bigdata);
                set(GUI.hcurrentfrequency, 'String', sprintf('%6.1f kHz', fr(j)));

                if (~strcmp(cmd, 'tone_test'))
                    set(GUI.hstimfilename, 'string', fnamep);
                else
                    set(GUI.hstimfilename, 'String', 'Tone Test (no file)');
                end
                
                for i = 1:length(spllist)
                    set(GUI.hcurrentSPL, 'String', sprintf('%3.1f dB', spllist(i)));
                    [DATA, STIM, HW, err] = tone_map(cmd, fr(j), spllist(i), HW, CALIBRATION, STIM, DATA, PLOTS, GUI); % main tone abr routine
                    if (err ~= 0)
                        err = 1;
                        return;
                    end
                    if ~strcmp(cmd, 'tone_test')
                        maxr(i, j) = max(abs(DATA.DATAa(1, :)')); % measure the signal
                        tmax = max(max(maxr));
                        maxry = 5 * maxr / tmax;
                        %                fprintf(2, 'maxr(i,j): %g   tmax: %g\n', maxr(i,j), tmax);
                        refreshdata(PLOTS.responsemap, 'caller');
                        drawnow
                    end

                    DATA.CDATAp = [DATA.CDATAp; DATA.DATAp]; % append the positive data
                    DATA.CDATAn = [DATA.CDATAn; DATA.DATAn]; % and the negative data
                    DATA.C_CHDATA = [DATA.C_CHDATA; DATA.CHDATA];
                    %   fprintf(1, 'fr = %6.1f, j = %d, spl = %6.1f  i = %d\n', fr(j), j, spls(i), i);
                    %               pause(0.1); % brief delay between frequencies.
                end

                GrandDatap{j} = DATA.CDATAp;
                GrandDatan{j} = DATA.CDATAn;
                DATA.CDATAp = DATA.CDATAp';
                DATA.CDATAn = DATA.CDATAn';
            %    DATA.C_CHDATA = DATA.C_CHDATA';

                if (~strcmp(cmd, 'tone_test'))
                    data_struct_p = DATA.CDATAp; % can't save a part of a structure?
                    data_struct_n = DATA.CDATAn;
                    save(fnamep, 'data_struct_p', '-ascii', '-tabs');
                    save(fnamen, 'data_struct_n', '-ascii', '-tabs');
                    % make a structure with all the data and parameters -
                    % even better!
                    bigdata = struct('CAL', CALIBRATION, 'DATA', DATA, 'STIM', STIM, 'HARDWARE', HW);
                    save(fname_bigdata, 'bigdata', '-mat');
                end

            end

            spl = spllist';

            if ~strcmp(cmd, 'tone_test')
                save(fnamei, 'spl', '-ascii', '-tabs'); % save intensity list also
                save(fnamef, 'fr', '-ascii', '-tabs'); % and frequency list
            end

            set(GUI.hstat, 'String', 'Done');
            last_data = cmd;


        otherwise
    end

end

%----------------------------------------------------------------------------
%****************************************************************************
%----------------------------------------------------------------------------

% one_click function to generate clicks and collect abr results This routine
% collects data for one stimulus condition (SPL) 9/23/03 P. Manis
%

function [DATA, STIM, err] = one_click(spl, mode, HW, STIM, CALIBRATION, DATA, PLOTS, GUI)
    % take a single parameter (spl, mode) click abr
    % requires call to getStimParams first.
    % if mode is 'test', we run in a reduced acquisition format
    %

    % fprintf(2, 'Calling one_click\n');
    err = 0;
    STIM.rate = 1.0 / STIM.sample_freq; % rate is in sec per point (recording side).
    STIM.click_dur = 0.1; % default is 100 microseconds (0.1)
    keep_reference = 1;
    DATA.REFERENCE = [];

    if strcmp(mode, 'test')
        oldsweeps = STIM.NSweeps;
        STIM.NSweeps = 10;
        oldsps = STIM.StimPerSweep;
        STIM.StimPerSweep = 40;
        oldipi = STIM.ipi;
        STIM.ipi = 25;
        set(GUI.hsr, 'String', num2str(STIM.NSweeps));
        STIM = updateStimParams(STIM, GUI);
    end

    cnp = STIM.StimPerSweep;

    if strcmp(HW.HARDWARE, 'NI')
        STIM.NIFreq = 500000; % express in sec per point
    else
        STIM.NIFreq = 44100;
    end

    STIM.delay = STIM.click_delay;
    STIM = updateStimParams(STIM, GUI);

    [STIM.wave, STIM.clock] = click(CALIBRATION.SPLCAL.click_amp, STIM.click_delay, STIM.click_dur, ...
        STIM.NIFreq, STIM.ipi, cnp, STIM.Alternate); % convert rate to usec per point
    attn = CALIBRATION.SPLCAL.maxclick - spl;

    if (attn < 0)
        fprintf(1, 'Requesting sound louder than available with this speaker\nSetting to 0 attn\n');
        attn = 0;
    end

    clear_plots(PLOTS, STIM, DATA, keep_reference);
    fl = get_SignalPlotFlag(1);
    STIM.Monitor = get_SignalPlotFlag(2);

    if fl == 1
        ts = (0:length(STIM.wave) - 1) * 1000 * STIM.clock;
        maxt = (length(STIM.wave) - 1) * 1000 * STIM.clock; % in msec
        %     figure; plot(ts,PLOTS.signal1);
        plot(PLOTS.signal1, ts, STIM.wave, 'color', 'blue');
        set(PLOTS.signal1, 'XLim', [0 maxt]);
        drawnow;
    end

    [data, STIM, chdata, err] = acquire4('attn', HW, STIM, PLOTS, GUI, CALIBRATION, attn);
    HW = set_attn(HW, -1);

    if (err == 0)
        [DATA, davg] = final_plot(data, STIM, PLOTS, DATA, keep_reference);
        DATA.DATAp = data(1, :);
        DATA.DATAn = data(2, :);
        DATA.DATAa = davg;
        DATA.CHDATA = chdata;
        pause(2.0);
    else
        fprintf(2, 'acquire4 err: %d\n', err);
        return
    end

    if keep_reference && err == 0
        DATA.REFERENCE = davg;
    end

    if strcmp(mode, 'test')
        STIM.NSweeps = oldsweeps;
        STIM.ipi = oldipi;
        STIM.StimPerSweep = oldsps;
        STIM = updateStimParams(STIM, GUI);
    end

end

%----------------------------------------------------------------------------
%****************************************************************************
%----------------------------------------------------------------------------

function [DATA, STIM, HW, err] = tone_map(mode, freq, spl, HW, CALIBRATION, STIM, DATA, PLOTS, GUI)

    err = 0;
    STIM.rate = 1 / STIM.sample_freq; % rate is in sec per point.

    if strcmp(HW.HARDWARE, 'NI')
        STIM.NIFreq = 500000; % express in sec per point
    else
        STIM.NIFreq = 44100;
    end

    if ~isempty(HW.AO)

        if (STIM.Alternate)
            tnp = STIM.StimPerSweep;
        else
            tnp = STIM.StimPerSweep;
        end

    else
        tnp = STIM.StimPerSweep;
    end

    if strcmp(mode, 'tone_test')
        DATA.REFERENCE = [];
        oldsweeps = STIM.NSweeps;
        STIM.NSweeps = 5;
        oldsps = STIM.StimPerSweep;
        STIM.StimPerSweep = 40;
        oldipi = STIM.ipi;
        STIM.ipi = 25;
        set(GUI.hsr, 'String', num2str(STIM.NSweeps));
    end

    STIM.delay = STIM.tone_delay;
    STIM = updateStimParams(STIM, GUI);
    % interpolate to get the attenuation at the requested frequency
    [splatF] = soundfuncs.spl_at_f(CALIBRATION.SPKR.Freqs, CALIBRATION.SPKR.maxdB, freq);

    % splatF=interp1(CALIBRATION.SPKR.Freqs, CALIBRATION.SPKR.maxdB, freq, 'spline');
    attn = splatF - spl;

    if (attn <= 0)
        fprintf(1, 'Requesting sound louder than available with this speaker\nSetting to 0 attn\n');
        return;
    end

    if strcmp(mode, 'tone_test')
        disp(splatF)
        disp(freq)
        disp(attn)
        disp(spl)
        disp(CALIBRATION.SPKR.CalAttn)
    end

    % fprintf(2, 'mode: %s', mode)
    if strcmp(mode, 'tone_info')
        disp(CALIBRATION.SPKR)
        fprintf(2, 'Freq: %8.2f  SPL: %8.2f  Attn: %6.1f\n', freq, spl, attn);
        err = 2;
        return;
    end

    switch (mode)
        case {'tone_abr', 'test', 'tone_test'}
            [STIM.wave, STIM.clock] = tonepip(CALIBRATION.SPLCAL.ToneVoltage, freq, STIM.tone_delay, ...
                5, 0.5, 0, STIM.NIFreq, STIM.ipi, tnp, STIM.Alternate); % convert rate to usec per point
        case {'tone_mapping', 'tone_info'}
            sel = get(GUI.hToneMapDuration, 'value');
            slist = str2num(get(GUI.hToneMapDuration, 'String')); %#ok<ST2NM>
            tdur = slist(sel);
            [STIM.wave, STIM.clock] = tonepip(CALIBRATION.SPLCAL.ToneVoltage, freq, STIM.tone_delay_mapping, tdur, 2.5, 0, STIM.NIFreq, STIM.ipi, STIM.StimPerSweep, 0); % convert rate to usec per point
        otherwise
            return;
    end

    clear_plots(PLOTS, STIM, DATA);
    fl = get_SignalPlotFlag(1);
    STIM.Monitor = get_SignalPlotFlag(2);

    if fl == 1
        ts = (0:length(STIM.wave) - 1) * 1000 * STIM.clock;
        maxt = (length(STIM.wave) - 1) * 1000 * STIM.clock; % in msec
        plot(PLOTS.signal1, ts, STIM.wave, 'color', 'blue');
        set(PLOTS.signal1, 'XLim', [0 maxt]);
        set(PLOTS.signal1, 'clipping', 'on');

        drawnow;
    end
    [data, STIM, chdata, err] = acquire4('attn', HW, STIM, PLOTS, GUI, CALIBRATION, attn);
    HW = set_attn(HW, -1);

    if (err ~= 0)
        return
    end

    if (err == 0)
        [DATA, davg] = final_plot(data, STIM, PLOTS, DATA);
        DATA.DATAp = data(1, :);
        DATA.DATAn = data(2, :);
        DATA.DATAa = davg;
        DATA.CHDATA = chdata;
    end

    if (strcmp(mode, 'tone_test'))
        STIM.NSweeps = oldsweeps;
        STIM.ipi = oldipi;
        STIM.StimPerSweep = oldsps;
        STIM = updateStimParams(STIM, GUI);
    end

end

function [DATA, davg] = final_plot(data, STIM, PLOTS, DATA, varargin)
    %tr = (0:length(data)-1)*STIM.rate*1000; % express rate in usec

    if (STIM.Alternate)
        davg = (data(1, :) + data(2, :)) / 2.0;
        dsum = (data(1, :) - data(2, :)) / 2.0;
    else
        davg = data(1, :); % already the average
        dsum = davg;
    end

    if nargin > 0
        DATA.REFERENCE = davg;
    end

    hold on;
    plot(PLOTS.data, STIM.ACQPars_tb, davg, 'k-');
    %md = max(abs(dsum));
    %trace_scale(md, PLOTS.data);
end

function [HW] = stop_acq(HW, status) 
    HW.RP.SoftTrg(0); % invoke(RP, 'softtrg', 0);
    HW.RP.Halt; % invoke(RP, 'halt');
    set_attn(HW, 120);
    %queueOutputData(AO, 0); % make sure output is really back to 0
    if strcmp(status, 'Stop')
        HW.IN_ACQ = false;
        set_status('Stopped');
    elseif strcmp(status, 'Pause')
        set_status('Paused')
    end
end


