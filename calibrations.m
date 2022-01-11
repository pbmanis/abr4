function calibrations( cmd, correctCal)
% This function handles the calibrations. There are two parts, depending on
% the command:
% 'Calibrate' records a range of frequencies at a fixed intensity
% and measures the microphone output, then calculates the SPL. The result is stored in a file
% 'frequencies_speakername.cal'. as a matlab structure.
% The other mode tests the microphone sensitivity (using a reference source setup)
% the results are stored in a file 'microphone_model#serialno.cal', for reference
% for future speaker calibrations.

% This code is part of ABR4, and references globals and GUI objects.
% 2010-2022 Paul B. Manis and lab...

global SPLCAL STIM SPKR STOP

STOP = 0;

if strcmp(cmd, 'calibrate')
    [Speaker, Mic] = getSpeakerMic(); % read from the gui
    calmode = 1;
    % Do not sample at 40 or 80 kHz - power supply noise interferes with measurement.
    switch (Speaker)
        case {'EC1', 'ES1'}
            spkr_freq = 1000*[1,2,3,4,5,6,7,8,9,10,12,14,16,20,22,24,26,28,30,32,36,38,42,44,48,52,56,60,64]; % speaker frequency list
        case 'MF1'
   %          spkr_freq = 1000*[87, 84, 82, 78, 76, 72, 68, 64, 60, 56, 52, ...
    %             48, 44, 42, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 11,...
    %             10,9, 8, 7, 6, 5, 4, 2, 1]; % speaker frequency list
%              spkr_freq = 1000*[87, 71.7, 58.8, 48, 37.5, 32.4, 28.0, 24.0, 20.5, 16, 14.0, 12.0, 10.0, 8, 6.0, 4, 3.0, 2, 1]; % speaker frequency list
%            spkr_freq = 1000*[ 12, 20.0, 35.0, 45.0, 80.0]; % speaker frequency list for quick testing.
%            spkr_freq = 1000*logspace(log10(0.5), log10(192.513/2.2), 27);
%              spkr_freq = [12500.0];
              spkr_freq = 1000*[20, 16, 12.5, 10, 8, 6.3, 5, 4, 3.15, 2.5, 2, 1.6, 1.25, 1.0];
             
        otherwise
            fprintf(2, 'Speaker type not known\n');
            return;
    end
else  % 'microphone' or 'microphone104'
    calmode = 0; % get microphone data
end

% hstat = findobj('tag', 'abr_status');
h_test = findobj('tag', 'TestCalc');
% h_mic = findobj('tag', 'ABR_Microphone');
set_attn(-1);
hf = findobj('tag', 'CalibrateFigure');
if isempty(hf)
    hf = figure;
    set(hf, 'tag', 'CalibrateFigure');
end
figure(hf);
clf;
checkSPL = 75.0; % level to test for correct calibration at all freqs
switch calmode
    case 1 % generate bank of frequencies to test
        err = verifyConfiguration('calibrate', cmd);
        if err == 1
            return;
        end
        fprintf(2, '\n==========================================\n');
        if correctCal == 0
            fprintf(2, 'Running Speaker Calibration on %s\n', date);
        else
            fprintf(2, 'Checking Speaker Calibration on %s\n', date);
        	CAL_L = load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file
            CAL = CAL_L.CAL;
            disp(CAL);
        end
        fprintf (2, 'Speaker: %s Microphone: %s Attenuation: %7.1f dB Reference level: %7.2f dBSPL\n', Speaker, Mic, SPKR.attn, checkSPL);
        MICL = load(sprintf('microphone_%s.cal', Mic), '-mat'); % get current microphone calibration file
        MIC = MICL.MIC;
        disp(MIC)
        amp = zeros(length(spkr_freq), 1);
%         theta = zeros(length(spkr_freq), 1);
%         M = zeros(length(spkr_freq), 1);
%         p_3a = zeros(length(spkr_freq), 1);
        dbspl = zeros(length(spkr_freq), 1);
        dbspl_bp = zeros(length(spkr_freq), 1);
        amnf = zeros(length(spkr_freq), 1);
        dbspl_nf = zeros(length(spkr_freq), 1);
        Vrms = zeros(length(spkr_freq), 1);
        Vrms_bp = zeros(length(spkr_freq), 1);
        Vrms_nf = zeros(length(spkr_freq), 1);
        attndB = zeros(length(spkr_freq), 1)+SPKR.attn;
        maxdB = zeros(length(spkr_freq), 1);
        recordDuration = 1.0; % seconds
        nRecordPoints = floor(STIM.sample_freq*recordDuration*1.5);
        ts1 = floor(STIM.sample_freq*0.1);
        ts2 = floor(STIM.sample_freq*0.9); % delay 100 msec, end at 900 msec
        
        fprintf(1, 'Recording Parameters: TraceDur: %7.1f s  points: %d STIM_sampleFreq: %9.3f Hz\n', ...
            recordDuration, nRecordPoints, STIM.sample_freq);
        rise_fall = 5.0; % msec
%         stopflag = 0;
%       design a notch filter at 40.0 kHz
        notchfilter = designfilt('bandstopiir', 'FilterOrder', 20, ...
            'HalfPowerFrequency1', 39000.0, 'HalfPowerFrequency2', 41000.0, ...
            'SampleRate', STIM.sample_freq);
%         fvtool(notchfilt);
%         return;
        trec = (0:1/STIM.sample_freq:(nRecordPoints-1)/STIM.sample_freq);
        for i = 1:length(spkr_freq)
            if check_stop(0) == 1
                fprintf(2, "Checkstop hit");
                return;
            end
            [STIM.wave, STIM.clock] = ...
                tonepip(SPLCAL.ToneVoltage, spkr_freq(i), ...
                0.0, recordDuration*1000., rise_fall, 0, STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
            if i == 1
                fprintf(1, 'Tone Parameters: TraceDur: %7.1f s  points: %d STIM.NIFreq: %9.3f Hz\n', ...
                    length(STIM.wave)/STIM.NIFreq, length(STIM.wave),  STIM.NIFreq);
                tstim = (0:1/STIM.NIFreq:(length(STIM.wave)-1)/STIM.NIFreq);
                fprintf(1, "Mic cal date: %s     Mic RefSig: %7.2f Mic mVPerPa:%7.2f\n", ...
                    MIC.Date, MIC.RefSig, MIC.mVPerPa);
                % column headers
                fprintf(2, '%8s\t%8s\t%7s  %6s  %6s %7s %7s %8s %7s %7s\n', ...
                    'F(Hz)', 'maxfreq', 'Mic(mV)','BP(mV)', 'NF(mv)', 'dBSPL','dbBP', ...
                    'Attn(dB)', 'Max dB', 'NF(dB)');
                
            end
            if correctCal == 1
                splatF=interp1(CAL.Freqs, CAL.dBSPL, spkr_freq(i), 'spline'); 
                attn = splatF - checkSPL + SPKR.attn;
                if(attn < 0)
                    attn = 0.0;
                end
            else
                attn = SPKR.attn;
            end
            attndB(i) = attn;
            set_attn(attn);
            [~, ch2, err] = calstim(nRecordPoints); % get the data...
            if err == 1
                fprintf(2, "Calibrations: calstim error");
                return;
            end
            ch2 = filter(notchfilter, ch2);
            set_attn(-1);
            [~, ch2nf, ~] = calstim(nRecordPoints); % make a noise floor measurement
            ch2nf = filter(notchfilter, ch2nf);


            bpfreqs = splfuncs.octave_calc(spkr_freq(i), 8, STIM.sample_freq);
%             fprintf(1, "spkr: %f  bp: %f  to %f\n", spkr_freq(i), bpfreqs(1), bpfreqs(2));
            ybp = bandpass(ch2, bpfreqs, STIM.sample_freq, ...
                'StopbandAttenuation', 60, "Steepness", 0.9);
            ynf = bandpass(ch2nf, bpfreqs, STIM.sample_freq,...
                'StopbandAttenuation', 60, "Steepness", 0.9);
%             fprintf(1, "bp calculated");
            Vrms_bp(i) = rms(ybp);
%             [~, imax_fft] = max(ybp);

            [amp_cosinor, fr_cosinor] = compute_cosinors(spkr_freq(i), trec(ts1:ts2), ch2(ts1:ts2));
            [~, k] = max(amp_cosinor);
            mfreq = fr_cosinor(k);
            amp(i) = amp_cosinor(k);

            [amnf(i), ~] = compute_cosinors([mfreq], trec, ch2nf); %#ok<NBRAK> 
%             fprintf(1, "cosinors calculated")
            [bp_data, bp_freqs] = periodogram(ybp,rectwin(length(ybp)),length(ybp),STIM.sample_freq);
            [nf_data, nf_freqs] = periodogram(ynf,rectwin(length(ynf)),length(ynf),STIM.sample_freq);
%             fprintf(1, "spectra calculated");
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
            
            % disp 'calc done'
            Vrms_bp(i) = Vrms_bp(i);
            Vrms(i) = amp(i)/sqrt(2); % convert cosinor to RMS
            Vrms_nf(i) = amnf(i); % rms(ynf); 
            % old calculation, based on cosinor amplitudes
%             dbspl(i) = MIC.RefSig + 20*log10(amp(i)/MIC.Vrms); %#ok<NODEF>
%             dbsplnf(i) = MIC.RefSig + 20*log10(2*amnf(i)/MIC.Vrms);
%           new calculation based on rms in narrow window around the
%           stimulus frequency.
%             fprintf(2, 'Vrms: %10.6e  Vref: %9.6f mic gain: %.1f\n', Vrms(i), MIC.Vref_bp, MIC.Gain);
            dbspl(i) =  splfuncs.compute_spl(Vrms(i), MIC);
            dbspl_bp(i) = splfuncs.compute_spl(Vrms_bp(i), MIC); 
            dbspl_nf(i) = splfuncs.compute_spl(Vrms_nf(i), MIC);
            maxdB(i) = dbspl(i) + attndB(i);
            fprintf(2, '%8.1f\t%8.1f\t%7.3f\t%7.3f\t%7.3f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n', ...
                spkr_freq(i), mfreq, 1000*Vrms(i), 1000*Vrms_bp(i), 1000*Vrms_nf(i), ...
                dbspl(i),dbspl_bp(i), attndB(i), maxdB(i), dbspl_nf(i));
        end
        spl_freqs = (min(spkr_freq):2000:max(spkr_freq));
        SPLatF = interp1(spkr_freq, dbspl, spl_freqs, 'spline'); 
        SPLatF_bp = interp1(spkr_freq, dbspl_bp, spl_freqs, 'spline'); 
        subplot(2,2,4);
        plot(spl_freqs, SPLatF, 'g-');
        hold on;
        plot(spl_freqs, SPLatF_bp, 'c-');
        stem(spkr_freq, dbspl, 'go-', 'fill');
        stem(spkr_freq, dbspl_bp, 'co-');
        hold on;
        sp = [spkr_freq fliplr(spkr_freq)];
        ndp = [dbspl_nf' 0*fliplr(dbspl_nf)'];
        ndp(ndp < 0) = 0; % clip at 0 dB spl..
        patch(sp, ndp, [0.6, 0.1, 0.1]);
        set(gca, 'XScale', 'log');
        if correctCal == 0 % then we should save the current measured values
            CAL.RefSPL = SPLCAL.maxtones;
            CAL.Freqs = spkr_freq;
            CAL.maxdB = maxdB;
            CAL.dBSPL = dbspl;
            CAL.dBSPL_bp = dbspl_bp;
            CAL.dBSPL_nf = dbspl_nf;
            CAL.Vmeas = Vrms;
            CAL.Vmeas_bp = Vrms_bp;
            CAL.Gain = 20; % dB setting microphone amplifier gain
            CAL.CalAttn = SPKR.attn; % attenuator setting at which calibration was done
            CAL.Speaker = Speaker;
            CAL.Microphone = Mic;
            CAL.Date = date;
            CAL.DateTime = datetime();
            calname = sprintf('frequency_%s_%s.cal', Speaker, date);
            save(calname, 'CAL');
            calname_short = sprintf('frequency_%s.cal', Speaker);
            save(calname_short, 'CAL');
            fprintf(2, 'Values saved in %s and %s\n', calname, calname_short);
        else
            CHKCAL.RefSPL = SPLCAL.maxtones;
            CHKCAL.Freqs = spkr_freq;
            CHKCAL.maxdB = maxdB;
            CHKCAL.dBSPL = dbspl;
            CHKCAL.dBSPL_nf = dbspl_nf;
            CHKCAL.dBSPL_bp = dbspl_bp;
            CHKCAL.Vmeas = Vrms;
            CAL.Vmeas_bp = Vrms_bp;
            CHKCAL.attndB = attndB;
            CHKCAL.Gain = 20; % dB setting microphone amplifier gain
            CHKCAL.CalAttn = SPKR.attn; % attenuator setting at which calibration was done
            CHKCAL.Speaker = Speaker;
            CHKCAL.Microphone = Mic;
            CHKCAL.Date = date;
            CHKCAL.DateTime = datetime();
            save(sprintf('chk75db_%s.cal', Speaker), 'CHKCAL');
            fprintf(2, 'Calibration Check saved in chk75db.cal\n');
        end
        
        fprintf(2, '==========================================\n\n');
        
    case 0 % microphone calibration at 1 kHz...
        err = verifyConfiguration('calibrate', cmd);
        if err == 1
            return;
        end
        testmode = get(h_test, 'Value');
        [Speaker, Mic] = getSpeakerMic();
        
        fprintf (1, '\n\nSpeaker: %s      Selected Microphone: %s \n', Speaker, Mic);
        testfreq = 5000.0;  % for daq channel, out and away from calibration source frequency.
        set_attn(120.0);  % no extraneous sound
        [STIM.wave, STIM.clock] = tonepip(SPLCAL.ToneVoltage, testfreq, ...
            0.0, 10, 0.0, 0, ...  % delay dur rf phase ...
            STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
        recordDuration = 1.0; % seconds
        nRecordPoints = floor(recordDuration*STIM.sample_freq);
        [~, ch2, ~] = calstim(nRecordPoints); % get the data...
        fprintf(1, "Recording completed\n");
        subplot(4,1,1);
        tstim = (0:1.0/STIM.NIFreq:(length(STIM.wave)-1.0)/STIM.NIFreq);
        plot(tstim', STIM.wave', 'b-');
        subplot(4,1,2);
        trec = (0:1/STIM.sample_freq:(nRecordPoints-1)/STIM.sample_freq);
        if testmode
            ch2 = 0.05*sin(trec*1000*2.0*pi);
            ch2 = ch2-mean(ch2);
        end
        plot(trec', ch2', 'r-');
        hold on;
        %  plot(trec', ch1', 'g-');
        nominal_freq = 1000.0;  % nominal calibrator frequency
        fr = (0.95*nominal_freq:1:1.05*nominal_freq); % not sure of exact output frequency, test several
        [amp_c, fr_c] = compute_cosinors(fr, trec, ch2);
        [~, imax_c] = max(amp_c);
        
        fprintf(1, "bandpass start... ");
        [ybp] = bandpass(ch2, [940, 1020], STIM.sample_freq, 'StopbandAttenuation', 60, "Steepness", 0.8);
        plot(trec, ybp, 'k-', 'LineWidth', 0.5);
        fprintf(1, '  bpdone ...');
        [bp_data, bp_freqs] = periodogram(ybp,rectwin(length(ybp)),length(ybp), STIM.sample_freq);
        fprintf(1, "bandpass and periodgram done");
        subplot(4,1,3);
        plot(bp_freqs, bp_data, 'r-');
        xlim([900, 1100]);
        [~, imax_bp] = max(bp_data);
        
%         pad = 1;
%         [amp_fft, fr_fft, x_pks, y_pks] = compute_fft(ch2, STIM.sample_freq, pad);
%         [fmin, ~] = min(fr);
%         [fmax, ~] = max(fr);
%         fk = find((fr_fft > fmin) & (fr_fft <= fmax));
%         subplot(2,1,2);
%         plot(fr_fft, amp_fft, 'b-');
        subplot(4,1,4);
        plot(fr_c, amp_c, 'm-');
%         plot(x_pks, y_pks, 'ro');
        xlim([900, 1100]);

%         [max_amp_fft, imax_fft] = max(amp_fft(fk));  % mostly for max freq position
        %         MIC.Vref = max_amp;
        MIC.Gain = 20.0; % dB setting
        if strcmp(cmd, 'microphone')
            MIC.RefSig = 94.0; % dB SPL for tone
        elseif strcmp(cmd, 'microphone104')
            MIC.RefSig = 104.0;
        end
        MIC.Vrms = rms(ch2);
        MIC.Vref_c = amp_c(imax_c)/sqrt(2);  % convert to RMS
        MIC.Vref_bp = rms(ybp);
        MIC.Microphone = Mic;
        MIC.Date = date;
        MIC.dBPerVPa = splfuncs.compute_dBPerVPa(MIC.Vref_bp, MIC);
        MIC.mVPerPa = splfuncs.compute_mVPerPa(MIC);
        
        fprintf(2, '\nMIC Vref (RMS): %12.6f Vrms(V) (raw trace)\n', MIC.Vrms);
        fprintf(2, 'MIC Vref (FFT): %12.6f Vmax(V) at %12.3f Hz (bandpassed), \n', MIC.Vref_bp, bp_freqs(imax_bp));
        fprintf(2, 'MIC Vref (cos): %12.6f Vmax(V) at %12.3f Hz, \n', MIC.Vref_c, fr_c(imax_c));
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
            save(sprintf('microphone_%s_%s.cal', Mic, date), 'MIC');
            save(sprintf('microphone_%s.cal', Mic), 'MIC');
        end
    otherwise
        fprintf(2, 'Calibration mode is not known.\n');
        
end
return;
end

function [amp, fr] = compute_cosinors(fr, trec, ch2)
    amp = zeros(length(fr), 1);
    theta = zeros(length(fr), 1);
    M = zeros(length(fr), 1);
    p_3a = zeros(length(fr), 1);
    for i = 1:length(fr)
        [amp(i), theta(i), M(i), p_3a(i)] = cosinor(trec', ch2', 2.0*pi*fr(i), 0.05);
    end

    return
end

function winplot(x, y, ts1, ts2, wincolor, outsidecolor)
plot(x(1:ts1), y(1:ts1), outsidecolor);
hold on
plot(x(ts1:ts2), y(ts1:ts2), wincolor);
plot(x(ts2:end), y(ts2:end), outsidecolor);
end

