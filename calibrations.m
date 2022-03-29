function [HW, STIM] = calibrations(cmd, check_calibration, HW, CALIBRATION, STIM, GUI)
% This function handles the sound system calibrations.
% There are two parts, depending on the argument passed in 'cmd':
% 1.  'calibrate' records a range of frequencies at a fixed intensity
% and measures the microphone output, then calculates the SPL. The result is stored in a file
% 'frequencies_speakername.cal'. as a matlab structure.
% The other mode tests the microphone sensitivity (using a reference source setup)
% the results are stored in a file 'microphone_model#serialno.cal', for reference
% for future speaker calibrations.

% This code is part of ABR4. It references GUI objects.
% 2010-2022 Paul B. Manis and lab...

if strcmp(cmd, 'calibrate')
    [Speaker, Mic] = getSpeakerMic(); % read from the gui
    calmode = 1;
    % Do not sample at 20, 40 or 80 kHz - power supply noise interferes with measurement.
    switch (Speaker)
        case {'EC1', 'ES1'}
            spkr_freq = 1000*[1,2,3,4,5,6,7,8,9,10,12,14,16,22,24,26,28,30,32,36,38,42,44,48,52,56,60,64]; % speaker frequency list
        case 'MF1'
           % Full:
           spkr_freq = 1000*[87, 84, 82, 78, 76, 72, 68, 64, 60, 56, 52, ...
                         48, 44, 42, 38, 36, 34, 32, 30, 28, 26, 24, 22, 18, 16, 14, 12, 11,...
                         10,9, 8, 7, 6, 5, 4, 3.2, 2.4, 2, 1.6, 1.2, 1]; % speaker frequency list
           % quick, broad
           %  spkr_freq = 1000*[87, 71.7, 58.8, 48, 37.5, 32.4, 28.0, 24.0, 20.5, 16, 14.0, 12.0, 10.0, 8, 6.0, 4, 3.0, 2, 1]; % speaker frequency list
           % Quick, high frequencies
           % spkr_freq = 1000*[ 12, 20.0, 35.0, 45.0, 80.0]; % speaker frequency list for quick testing.
           %
           % log spaced frequencies (too far separated for low frequencies)
           %            spkr_freq = 1000*logspace(log10(0.5), log10(192.513/2.2), 27);
          %
          % very quick, wide range
          %  spkr_freq = 1000*[ 72, 60, 48, 32, 16, 12, 8, 4, 2.5, 1 ];
          % quick, low frequencies.
          % spkr_freq = 1000*[20, 16, 12.5, 10, 8, 6.3, 5, 4, 3.2, 2.4, 2, 1.6, 1.2, 1.0];
            
        otherwise
            fprintf(2, 'Speaker type not known\n');
            return;
    end
elseif strcmp(cmd, 'microphone') || strcmp(cmd, 'microphone104')
    calmode = 0; % get microphone data
else
    fprintf(2, 'Call to calibrations with command "%s" is not recognized', cmd)
    return
    
end

% hstat = findobj('tag', 'abr_status');
h_test = findobj('tag', 'TestCalc');
% h_mic = findobj('tag', 'ABR_Microphone');
HW = set_attn(HW, -1);
hf = findobj('tag', 'CalibrateFigure');
if isempty(hf)
    hf = figure;
    set(hf, 'tag', 'CalibrateFigure');
end
figure(hf);
clf;
checkSPL = 75.0; % level to test for correct calibration at all freqs
switch calmode
    case 0 
        %------------------------------------------------------------------
        % perform microphone calibration using a sound source, nominally at
        % 1 kHz...
        %------------------------------------------------------------------
        err = verifyConfiguration('calibrate', cmd);
        if err == 1
            return;
        end
        testmode = get(h_test, 'Value');  % get the testing mode (calculations only or full sound calibration)
        [Speaker, Mic] = getSpeakerMic();
        fprintf (1, '\n\nSpeaker: %s      Selected Microphone: %s \n', Speaker, Mic);

        testfreq = 5000.0;  % for daq channel, out and away from calibration source frequency.
        HW = set_attn(HW, 120.0);  % no extraneous sound
        [STIM.wave, STIM.clock] = tonepip(CALIBRATION.SPLCAL.ToneVoltage, testfreq, ...
            0.0, 10, 0.0, 0, ...  % delay dur rf phase ...
            STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
        recordDuration = 1.0; % seconds
        nRecordPoints = floor(recordDuration*STIM.sample_freq);
        [~, ch2, HW, ~] = calstim(nRecordPoints, HW, STIM); % get the data...
        fprintf(1, "   ... calibration sound recorded\n");
        % display the stimulus waveform
        subplot(4,1,1);
        tstim = (0:1.0/STIM.NIFreq:(length(STIM.wave)-1.0)/STIM.NIFreq);
        plot(tstim', STIM.wave', 'b-');
        % display the recorded waveform, with parts identified
        subplot(4,1,2);
        trec = (0:1/STIM.sample_freq:(nRecordPoints-1)/STIM.sample_freq);
        if testmode
            ch2 = 0.05*sin(trec*1000*2.0*pi);
            ch2 = ch2-mean(ch2);
        end
        plot(trec, ch2, 'r-');
        hold on;
        nominal_freq = 1000.0;  % nominal calibrator frequency
        fr = (0.95*nominal_freq:1:1.05*nominal_freq); % not sure of exact output frequency, test several
        [amp_cs, fr_cs] = compute_cosinors(fr, trec, ch2);
        [~, imax_cs] = max(amp_cs);
        
        fprintf(1, "Starting bandpass filtering calculation ... ");
        [ybp] = bandpass(ch2, [940, 1020], STIM.sample_freq, 'StopbandAttenuation', 60, "Steepness", 0.8);
        [bp_data, bp_freqs] = periodogram(ybp,rectwin(length(ybp)),length(ybp), STIM.sample_freq);
        
        plot(trec, ybp, 'k-', 'LineWidth', 0.5);
        fprintf(1, " Bandpass and periodgram done\n");
        subplot(4,1,3);
        plot(bp_freqs, bp_data, 'r-');  % bandpass spectrum.
        xlim([900, 1100]);
        [~, imax_bp] = max(bp_data); 
        subplot(4,1,4);
        plot(fr_cs, amp_cs, 'm-');  % cosinor calculated spectrum
        xlim([900, 1100]);
        MIC.Gain = 20.0; % dB setting
        if strcmp(cmd, 'microphone')
            MIC.RefSig = 94.0; % dB SPL for tone
        elseif strcmp(cmd, 'microphone104')
            MIC.RefSig = 104.0;
        end
        %
        % Microphone voltage extraction and calculations of transfer
        % factors
        MIC.Vrms = rms(ch2);
        MIC.Vref_c = amp_cs(imax_cs)/sqrt(2);  % convert to RMS
        MIC.Vref_bp = rms(ybp);
        MIC.Microphone = Mic;
        MIC.Date = date;
        MIC.dBPerVPa = soundfuncs.compute_dBPerVPa(MIC.Vref_bp, MIC);
        MIC.mVPerPa = soundfuncs.compute_mVPerPa(MIC);
        
        fprintf(2, '\nMIC Vref (RMS): %12.6f Vrms(V) (raw trace)\n', MIC.Vrms);
        fprintf(2, 'MIC Vref (BP) : %12.6f Vmax(V) at %12.3f Hz (bandpassed), \n', MIC.Vref_bp, bp_freqs(imax_bp));
        fprintf(2, 'MIC Vref (cos): %12.6f Vmax(V) at %12.3f Hz, \n', MIC.Vref_c, fr_cs(imax_cs));
        fprintf(1, '---------------------------------------------------------\n');
        switch Mic
            case '7012#39279'
                fprintf(1, 'Calibration information from Manufacturer, June 5, 2008\n');
                fprintf(1, '1/2\" 7012 mic (SN 39279):  %7.2f mV/Pa, -35.8 dB re 1V/Pa\n', 16.22);
            case '7016#9945'
                fprintf(1, 'Calibration information from Manufacturer, June 3, 2010\n');
                fprintf(1, '1/4\" 7016 mic (SN  9945):  %7.2f mV/Pa, -49.0 dB re 1V/Pa\n', 3.55);
            case '7016#10252'
                fprintf(1, 'Calibration information from Manufacturer, 2018 (7016 #10252) \n');
                fprintf(1, '1/4\" 7016 mic (SN  10252):  %7.2f mV/Pa, -48.3 dB re 1V/Pa\n', 3.85);
        end
        fprintf(2, 'Measured Transfer factor:   %7.2f mV/Pa, %5.1f dB re 1V/Pa\n', MIC.mVPerPa, MIC.dBPerVPa);
        fprintf(2, 'Measured with Mic Amp Gain = 20.0dB, Standard is %.1f dB, bandpassed\n', MIC.RefSig);
        fprintf(1, '--------------------------------------------------------\n');
        
        if ~testmode
            calfile = sprintf('calibration_history\\microphone_%s_%s.cal', Mic, date);
            save(calfile, 'MIC');
            save(sprintf('microphone_%s.cal', Mic), 'MIC');
        end
    

    case 1 % 
        %------------------------------------------------------------------
        % Perform a speaker calibration using a range of frequencies
        %------------------------------------------------------------------
        err = verifyConfiguration('calibrate', cmd);
        if err == 1
            return;
        end
        fprintf(2, '\n==========================================\n');
        if check_calibration == 0
            fprintf(2, 'Running Speaker Calibration on %s\n', date);
        else
            fprintf(2, 'Checking Speaker Calibration on %s\n', date);
            CAL_L = load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file
            CAL = CAL_L.CAL;
            disp(CAL);
        end
        fprintf (2, 'Speaker: %s Microphone: %s Attenuation: %7.1f dB Reference level: %7.2f dBSPL\n', ...
            Speaker, Mic, CALIBRATION.SPKR.CalAttn, checkSPL);
        MICL = load(sprintf('microphone_%s.cal', Mic), '-mat'); % get the current microphone calibration file
        MIC = MICL.MIC;

        amp_cs = zeros(length(spkr_freq), 1);  % cosinor amplitude raw
        amp_nf = zeros(length(spkr_freq), 1);  % noise floor amplitude raw
 
        dbspl_cs = zeros(length(spkr_freq), 1);  % cosinor db
        dbspl_bp = zeros(length(spkr_freq), 1);  % bandpassed db
        dbspl_nf = zeros(length(spkr_freq), 1);  % noise floor db
        
        Vrms_cs = zeros(length(spkr_freq), 1);  % cosinor rms measure mic voltage
        Vrms_bp = zeros(length(spkr_freq), 1); % bandpassed rms mic voltage
        Vrms_nf = zeros(length(spkr_freq), 1);  % noise floor rms mic voltage
        
        attndB = zeros(length(spkr_freq), 1)+CALIBRATION.SPKR.CalAttn;
        maxdB = zeros(length(spkr_freq), 1);
        
        if check_calibration == 1
            recordDuration = 0.5; % seconds
        else
            recordDuration = 0.2; % seconds
        end
        nRecordPoints = floor(STIM.sample_freq*recordDuration*1.5);
        ts1 = floor(recordDuration*STIM.sample_freq*0.1);
        ts2 = floor(recordDuration*STIM.sample_freq*0.9); % delay 100 msec, end at 900 msec
        fprintf(1, 'Recording Parameters: TraceDur: %7.1f s  points: %d STIM_sampleFreq: %9.3f Hz\n', ...
            recordDuration, nRecordPoints, STIM.sample_freq);
        rise_fall = 5.0; % msec

        % design a notch filter at 40.0 kHz
        notchfilter1 = designfilt('bandstopiir', 'FilterOrder', 20, ...
            'HalfPowerFrequency1', 39000.0, 'HalfPowerFrequency2', 41000.0, ...
            'SampleRate', STIM.sample_freq);
        % design a notch filter at 20.0 kHz
        notchfilter2 = designfilt('bandstopiir', 'FilterOrder', 20, ...
            'HalfPowerFrequency1', 19500.0, 'HalfPowerFrequency2', 20500.0, ...
            'SampleRate', STIM.sample_freq);
        %         fvtool(notchfilt);  % to view the bandpass of the filter we just made
        %         return;
        trec = (0:1/STIM.sample_freq:(nRecordPoints-1)/STIM.sample_freq);
        for i = 1:length(spkr_freq)
%             state = check_status(GUI);
%             if state == 'Stopped'
%                 fprintf(2, "Abort hit, terminating speaker calibration without saving data\n");
%                 return;
%             end
            [STIM.wave, STIM.clock] = ...
                tonepip(CALIBRATION.SPLCAL.ToneVoltage, spkr_freq(i), ...
                0.0, recordDuration*1000., rise_fall, 0, STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
            if i == 1  % print information and a header for table
                fprintf(1, 'Tone Parameters: TraceDur: %7.1f s  points: %d STIM.NIFreq: %9.3f Hz\n', ...
                    length(STIM.wave)/STIM.NIFreq, length(STIM.wave),  STIM.NIFreq);
                tstim = (0:1/STIM.NIFreq:(length(STIM.wave)-1)/STIM.NIFreq);
                fprintf(1, "Mic cal date: %s     Mic RefSig: %7.2f Mic mVPerPa:%7.2f\n", ...
                    MIC.Date, MIC.RefSig, MIC.mVPerPa);

                fprintf(2, '%8s\t%8s\t%7s  %6s  %6s %7s %7s %8s %7s %7s\n', ...
                    'F(Hz)', 'maxfreq', 'Mic(mV)','BP(mV)', 'NF(mv)', 'dBWB','dbBP', ...
                    'Attn(dB)', 'Max dB', 'dbNF');
            end
            if check_calibration == 1  % 
                %splatF=interp1(CAL.Freqs, CAL.dBSPL, spkr_freq(i), 'spline');
                splatF = soundfuncs.spl_at_f(CAL.Freqs, CAL.dBSPL_bp, spkr_freq(i)); % linear on log scale.
                attn = splatF - checkSPL + CALIBRATION.SPKR.CalAttn;
                if(attn < 0)
                    attn = 0.0;
                end
            else
                attn = CALIBRATION.SPKR.CalAttn;
            end
            attndB(i) = attn;
            HW = set_attn(HW, attn);
            [~, ch2, HW, err] = calstim(nRecordPoints, HW, STIM); % get the data...
            if err == 1
                fprintf(2, "Calibrations: calstim error");
                return;
            end
            ch2 = filter(notchfilter1, ch2);
            ch2 = filter(notchfilter2, ch2);
          
            HW = set_attn(HW, -1);
            [~, ch2nf, HW, ~] = calstim(nRecordPoints, HW, STIM); % make a noise floor measurement
            ch2nf = filter(notchfilter1, ch2nf);
            ch2nf = filter(notchfilter2, ch2nf);
            
            bpfreqs = soundfuncs.octave_calc(spkr_freq(i), 5, STIM.sample_freq);
            %             fprintf(1, "spkr: %f  bp: %f  to %f\n", spkr_freq(i), bpfreqs(1), bpfreqs(2));
            ybp = bandpass(ch2, bpfreqs, STIM.sample_freq, ...
                'StopbandAttenuation', 60, "Steepness", 0.9);
            ynf = bandpass(ch2nf, bpfreqs, STIM.sample_freq,...
                'StopbandAttenuation', 60, "Steepness", 0.9);
            %             fprintf(1, "bp calculated");
            Vrms_bp(i) = rms(ybp);
            Vrms_nf(i) = rms(ynf);
            
            [amp_cosinor, fr_cosinor] = compute_cosinors(spkr_freq(i), trec(ts1:ts2), ch2(ts1:ts2));
            [~, k] = max(amp_cosinor);
            mfreq = fr_cosinor(k);
            amp_cs(i) = amp_cosinor(k);
            Vrms_cs(i) = amp_cs(i)/sqrt(2); % convert cosinor to RMS
            
%             [amp_nf(i), ~] = compute_cosinors([mfreq], trec, ch2nf); %#ok<NBRAK>
%             Vrms_nf(i) = amp_nf(i)/sqrt(2);
 
            [bp_data, bp_freqs] = periodogram(ybp,rectwin(length(ybp)),length(ybp),STIM.sample_freq);
            [nf_data, nf_freqs] = periodogram(ynf,rectwin(length(ynf)),length(ynf),STIM.sample_freq);
 
            figure(hf);
            clf;
            % plot region of raw stimulus used
            subplot(2,2,1);
            plot(tstim, STIM.wave);
            % data in raw recording, same window
            subplot(2,2,2);
            winplot(trec, ch2, ts1, ts2, 'k-', 'b-');
            % bandpassed signals
            winplot(trec, ybp, ts1, ts2, 'r-', 'm-');
            winplot(trec, ynf, ts1, ts2, 'c-', 'y-');
            % spectra in the region of the signal
            subplot(2,2,3);
            plot(bp_freqs, bp_data, 'b-');
            hold on;
            plot(nf_freqs, nf_data, 'm-');
            xlim([spkr_freq(i)-1000.0, spkr_freq(i) + 1000.0])
            
         %   Vrms_bp(i) = Vrms_bp(i);
            % old calculation (2011-Jan 2022), based on cosinor amplitudes
            %             dbspl(i) = MIC.RefSig + 20*log10(amp(i)/MIC.Vrms); %#ok<NODEF>
            %             dbsplnf(i) = MIC.RefSig + 20*log10(2*amnf(i)/MIC.Vrms);
            %           new calculation based on rms in narrow window around the
            %           stimulus frequency.
            %             fprintf(2, 'Vrms: %10.6e  Vref: %9.6f mic gain: %.1f\n', Vrms(i), MIC.Vref_bp, MIC.Gain);
            dbspl_cs(i) = soundfuncs.compute_spl(Vrms_cs(i), MIC);
            dbspl_bp(i) = soundfuncs.compute_spl(Vrms_bp(i), MIC);
            dbspl_nf(i) = soundfuncs.compute_spl(Vrms_nf(i), MIC);
            maxdB(i) = dbspl_bp(i) + attndB(i);
            fprintf(2, '%8.1f\t%8.1f\t%7.3f\t%7.3f\t%7.3f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n', ...
                spkr_freq(i), mfreq, 1000*Vrms_cs(i), 1000*Vrms_bp(i), 1000*Vrms_nf(i), ...
                dbspl_cs(i),dbspl_bp(i), attndB(i), maxdB(i), dbspl_nf(i));
        end
        spl_freqs = (min(spkr_freq):2000:max(spkr_freq));
%         SPLatF = interp1(spkr_freq, dbspl, spl_freqs, 'spline');
%         SPLatF_bp = interp1(spkr_freq, dbspl_bp, spl_freqs, 'spline');
        SPLatF_cs = soundfuncs.spl_at_f(spkr_freq, dbspl_cs, spl_freqs);
        SPLatF_bp = soundfuncs.spl_at_f(spkr_freq, dbspl_bp, spl_freqs);
        subplot(2,2,4);
        plot(spl_freqs, SPLatF_bp, 'g-');
        hold on;
%         plot(spl_freqs, SPLatF_cs, 'c-');
%         stem(spkr_freq, dbspl_cs, 'co-', 'fill');
        stem(spkr_freq, dbspl_bp, 'go-');
        hold on;
        sp = [spkr_freq fliplr(spkr_freq)];
        ndp = [dbspl_nf' 0*fliplr(dbspl_nf)'];
        ndp(ndp < 0) = 0; % clip at 0 dB spl..
        patch(sp, ndp, [0.6, 0.1, 0.1]);
        set(gca, 'XScale', 'log');
        
        if check_calibration == 0 % Save the current measured values for reference
            CAL.RefSPL = CALIBRATION.SPLCAL.maxtones;
            CAL.Freqs = spkr_freq;
            CAL.maxdB = maxdB;
            CAL.dBSPL = dbspl_cs;
            CAL.dBSPL_bp = dbspl_bp;
            CAL.dBSPL_nf = dbspl_nf;
            CAL.Vmeas = Vrms_cs;
            CAL.Vmeas_bp = Vrms_bp;
            CAL.Gain = 20; % dB setting microphone amplifier gain
            CAL.CalAttn = CALIBRATION.SPKR.CalAttn; % attenuator setting at which calibration was done
            CAL.Speaker = Speaker;
            CAL.Microphone = Mic;
            CAL.Date = date;
            CAL.DateTime = datetime();
            calname = sprintf('calibration_history\\frequency_%s_%s.cal', Speaker, date);
            save(calname, 'CAL');
            calname_short = sprintf('frequency_%s.cal', Speaker);
            save(calname_short, 'CAL');
            fprintf(2, 'Speaker response saved in %s and %s\n', calname, calname_short);
        else
            CHKCAL.RefSPL = CALIBRATION.SPLCAL.maxtones;
            CHKCAL.Freqs = spkr_freq;
            CHKCAL.maxdB = maxdB;
            CHKCAL.dBSPL = checkSPL;
            CHKCAL.dBSPL_nf = dbspl_nf;
            CHKCAL.dBSPL_bp = dbspl_bp;
            CHKCAL.Vmeas = Vrms_cs;
            CHKCAL.Vmeas_bp = Vrms_bp;
            CHKCAL.attndB = attndB;
            CHKCAL.Gain = 20; % dB setting microphone amplifier gain
            CHKCAL.CalAttn = CALIBRATION.SPKR.CalAttn; % attenuator setting at which calibration was done
            CHKCAL.Speaker = Speaker;
            CHKCAL.Microphone = Mic;
            CHKCAL.Date = date;
            CHKCAL.DateTime = datetime();
            save(sprintf('chk75db_%s.cal', Speaker), 'CHKCAL');
            save(sprintf('calibration_history\\chk75db_%s_%s.cal', Speaker, date), 'CHKCAL');
            fprintf(2, 'Calibration Check saved in chk75db.cal (and in calibration_history)\n');
        end
        
        fprintf(2, '==========================================\n\n');
        
    otherwise
        fprintf(2, 'Calibration mode is not known.\n');
        
end
return;
end

function [amp, fr] = compute_cosinors(fr, trec, wave)
% Compute the the amplitudes at specific frequencies
% fr : frequencies to compute
% trec : time base for the waveform
% wave : the waveform
% returns: lists of amplitudes and frequences.
amp = zeros(length(fr), 1);
theta = zeros(length(fr), 1);
M = zeros(length(fr), 1);
p_3a = zeros(length(fr), 1);
for i = 1:length(fr)
    [amp(i), theta(i), M(i), p_3a(i)] = cosinor(trec', wave', 2.0*pi*fr(i), 0.05);
end
end

function winplot(x, y, ts1, ts2, wincolor, outsidecolor)
% Plot x,y within a time window (points) in one color, and
% the datapoinbs outside that window in another color.
% Uses the current plot
plot(x(1:ts1), y(1:ts1), outsidecolor);
hold on;
plot(x(ts1:ts2), y(ts1:ts2), wincolor);
plot(x(ts2:end), y(ts2:end), outsidecolor);
end

