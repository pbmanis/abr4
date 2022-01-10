function [ output_args ] = calibrations( cmd, correctCal)
% This funciton handles the calibrations. There are two parts, depending on
% the command:
% 'Calibrate' records a range of frequencies at a fixed intensity
% and measures the microphone output, then calculates the SPL. The result is stored in a file
% 'frequencies_speakername.cal'. as a matlab structure.
% The other mode tests the microphone sensitivity (using a reference source setup)
% the results are stored in a file 'microphone_model#serialno.cal', for reference
% for future speaker calibrations.

% This code is part of ABR4, and references globals and GUI objects.
% 2010-2022 Paul B. Manis and lab...

global SPLCAL STIM SPKR

if strcmp(cmd, 'calibrate')
    [Speaker, Mic] = getSpeakerMic(); % read from the gui
    calmode = 1;
    % Do not sample at 40 or 80 kHz - power supply noise interferes with measurement.
    switch (Speaker)
        case {'EC1', 'ES1'}
            spkr_freq = 1000*[1,2,3,4,5,6,7,8,9,10,12,14,16,20,22,24,26,28,30,32,36,38,42,44,48,52,56,60,64]; % speaker frequency list
        case 'MF1'
            spkr_freq = 1000*[100, 96, 92, 88, 84, 82, 78, 76, 72, 68, 64, 56, 48, 44, 42, 38, 36, 32, 28, 24, 20, 18, 16, 12, 8, 6, 4, 2, 1]; % speaker frequency list
          %  spkr_freq = 1000*[0.5, 1, 2, 4, 8, 12]; % speaker frequency list for quick testing.
        otherwise
            fprintf(2, 'Speaker type not known\n');
            return;
    end
else  % 'microphone' or 'microphone104'
    calmode = 0; % get microphone data
end

hstat = findobj('tag', 'abr_status');
h_test = findobj('tag', 'TestCalc');
h_mic = findobj('tag', 'ABR_Microphone');
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
        end
        fprintf (2, 'Speaker: %s Microphone: %s Attenuation: %7.1f dB Reference level: %7.2f dBSPL\n', Speaker, Mic, SPKR.attn, checkSPL);
        load(sprintf('microphone_%s.cal', Mic), '-mat'); % get current microphone calibration file
        amp = zeros(length(spkr_freq), 1);
        theta = zeros(length(spkr_freq), 1);
        M = zeros(length(spkr_freq), 1);
        p_3a = zeros(length(spkr_freq), 1);
        dbspl = zeros(length(spkr_freq), 1);
        amnf = zeros(length(spkr_freq), 1);
        dbsplnf = zeros(length(spkr_freq), 1);
        Vrms = zeros(length(spkr_freq), 1);
        attndB = zeros(length(spkr_freq), 1)+SPKR.attn;
        TraceDuration = 2.0; % seconds
        nRecordPoints = floor(STIM.sample_freq*TraceDuration*1.2);
        fprintf(1, 'Recording Parameters: TraceDur: %7.1f s  points: %d STIM_sampleFreq: %9.3f Hz\n', ...
            TraceDuration, nRecordPoints, STIM.sample_freq);
        rise_fall = 5.0; % msec
        %tonepip(SPLCAL.ToneVoltage, spkr_freq(1), ...
        %        0.0, 1000, rise_fall, 0, STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
        stopflag = 0;
        trec = (0:1/STIM.sample_freq:(nRecordPoints-1)/STIM.sample_freq);
        for i = 1:length(spkr_freq)
            if check_stop(0) == 1
                return;
            end
            [STIM.wave, STIM.clock] = tonepip(SPLCAL.ToneVoltage, spkr_freq(i), ...
                0.0, TraceDuration*1000., rise_fall, 0, STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
            if i == 1
                fprintf(1, 'Tone Parameters: TraceDur: %7.1f s  points: %d STIM.NIFreq: %9.3f Hz\n', ...
                    length(STIM.wave)/STIM.NIFreq, length(STIM.wave),  STIM.NIFreq);
                tstim = (0:1/STIM.NIFreq:(length(STIM.wave)-1)/STIM.NIFreq);
                fprintf(2, '%8s\t %8s\t %12s\t%12s\t%12s\t%12s\n', ...
                    'F(Hz)', 'maxfreq', 'Mic(mV)','dBSPL','Attn(dB)', 'Vrms(mV)'); %, spkr_freq(i), mfreq, amp(i), dbspl(i), attndB(i));
                
            end
            if correctCal == 1
                load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file
                splatF=interp1(CAL.Freqs, CAL.maxdB, spkr_freq(i), 'spline'); %#ok<NODEF>
                attn = splatF - checkSPL + SPKR.attn;
                if(attn < 0)
                    attn = 0.0;
                end
            else
                attn = SPKR.attn;
            end
            attndB(i) = attn;
            set_attn(attn);
            [ch1, ch2, err] = calstim(nRecordPoints); % get the data...
            if err == 1
                return;
            end
            set_attn(-1);
            [~, ch2nf, err] = calstim(nRecordPoints); % make a noise floor measurement
            figure(hf);
            clf;
            subplot(2,2,1);
            plot(tstim, STIM.wave);
            subplot(2,2,2);
            plot(trec, ch2, 'r-');
            hold on;
            %plot(trec, ch1, 'g-');
            ts1 = floor(0.15*length(trec));
            ts2 = floor(0.75*length(trec));
            plot(trec(ts1:ts2), ...
                ch2(ts1:ts2), 'g-');
            range = 1.0;
            % disp 'calculate'
            clippts = 2*floor((rise_fall/1000.)*STIM.sample_freq);
            fr = [spkr_freq(i)*range];
            [amp_cosinor, fr_cosinor] = compute_cosinors(fr, trec, ch2);
            subplot(2,2,3);
            periodogram(ch2(ts1:ts2),rectwin(length(ch2(ts1:ts2))),length(ch2(ts1:ts2)),STIM.sample_freq)
            [~, k] = max(amp_cosinor);
            mfreq = fr_cosinor(k);
            amp(i) = amp_cosinor(k);
            [amnf(i), ~] = compute_cosinors([mfreq], trec, ch2nf);
            
            % disp 'calc done'
            Vrms(i) = rms(ch2);
            dbspl(i) = MIC.RefSig + 20*log10(amp(i)/MIC.Vrms); %#ok<NODEF>
            dbsplnf(i) = MIC.RefSig + 20*log10(2*amnf(i)/MIC.Vrms);
            fprintf(2, '%8.2f\t%8.2f\t%12.3f\t%12.1f\t%12.1f\t%12.3f\n', ...
                spkr_freq(i), mfreq, 1000.*amp(i), dbspl(i), attndB(i), 1000.0*Vrms(i));
        end
        %clf;
        subplot(2,2,4);
        stem(spkr_freq, dbspl, 'ko-', 'fill');
        hold on;
        sp = [spkr_freq fliplr(spkr_freq)];
        ndp = [dbsplnf' 0*fliplr(dbsplnf)'];
        ndp(ndp < 0) = 0; % clip at 0 dB spl..
        patch(sp, ndp, [0.6, 0.1, 0.1]);
        set(gca, 'XScale', 'log');
        if correctCal == 0 % then we should save the current measured values
            CAL.RefSPL = SPLCAL.maxtones;
            CAL.Freqs = spkr_freq;
            CAL.maxdB = dbspl;
            CAL.Vmeas = 2 * amp;
            CAL.Gain = 20; % dB setting microphone amplifier gain
            CAL.CalAttn = SPKR.attn; % attenuator setting at which calibration was done
            CAL.Speaker = Speaker;
            CAL.Microphone = Mic;
            CAL.Date = date;
            calname = sprintf('frequency_%s_%s.cal', Speaker, date);
            save(calname, 'CAL');
            calname_short = sprintf('frequency_%s.cal', Speaker);
            save(calname_short, 'CAL');
            fprintf(2, 'Values saved in %s and %s\n', calname, calname_short);
        else
            CHKCAL.RefSPL = SPLCAL.maxtones;
            CHKCAL.Freqs = spkr_freq;
            CHKCAL.maxdB = dbspl;
            CHKCAL.Vmeas = 2 * amp;
            CHKCAL.attndB = attndB;
            CHKCAL.Gain = 20; % dB setting microphone amplifier gain
            CHKCAL.CalAttn = SPKR.attn; % attenuator setting at which calibration was done
            CHKCAL.Speaker = Speaker;
            CHKCAL.Microphone = Mic;
            CHKCAL.Date = date;
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
        [STIM.wave, STIM.clock] = tonepip(SPLCAL.ToneVoltage, testfreq, 0.0, 10, 0.0, 0, STIM.NIFreq, 10, 1, 0); % convert rate to usec per point
        nRecordPoints = floor(5.0*STIM.sample_freq);  % take 5 seconds of data
        [ch1, ch2, ~] = calstim(nRecordPoints); % get the data...
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
        [max_amp_c, imax_c] = max(amp_c);
        pad = 1;
        [amp_fft, fr_fft, x_pks, y_pks] = compute_fft(ch2, STIM.sample_freq, pad);
        [fmin, ~] = min(fr);
        [fmax, ~] = max(fr);
        fk = find((fr_fft > fmin) & (fr_fft <= fmax));
        subplot(2,1,2);
        plot(fr_fft, amp_fft, 'b-');
        hold on;
        plot(fr_c, amp_c, 'm-');
%         plot(x_pks, y_pks, 'ro');
        xlim([900, 1100]);

        [max_amp_fft, imax_fft] = max(amp_fft(fk));  % mostly for max freq position
        %         MIC.Vref = max_amp;
        MIC.Gain = 20.0; % dB setting
        if strcmp(cmd, 'microphone')
            MIC.RefSig = 94.0; % dB SPL for tone
        elseif strcmp(cmd, 'microphone104')
            MIC.RefSig = 104.0;
        end
        MIC.Vrms = rms(ch2);
        MIC.Vref_c = amp_c(imax_c);
        MIC.Vref_fft = max_amp_fft;
        MIC.Microphone = Mic;
        MIC.Date = date;
        MIC.dBPerVPa = 20.0*log10(MIC.Vrms) - MIC.Gain - (MIC.RefSig - 94.0); % sensitivity.
        MIC.mVPerPa = 1000*10.0^(MIC.dBPerVPa/20.0);
        
        fprintf(2, '\nMIC Vref (RMS): %12.6f Vrms\n', MIC.Vrms);
        fprintf(2, '\nMIC Vref (FFT): %12.6f Vmax at %12.3f Hz, \n', MIC.Vref_fft, fr_fft(imax_fft));
        fprintf(2, '\nMIC Vref (cos): %12.6f Vmax at %12.3f Hz, \n', MIC.Vref_c, fr_c(imax_c));
        fprintf(1, '---------------------------------------------------------\n');
        switch Mic
            case '7012#39279'
                fprintf(1, 'Calibration information from Manufacturer, June 5, 2008\n');
                fprintf(1, '1/2\" 7012 mic (SN 39279):  16.22 mV/Pa, -35.8 dB re 1V/Pa\n');
            case '7016#9945'
                fprintf(1, 'Calibration information from Manufacturer, June 3, 2010\n');
                fprintf(1, '1/4\" 7016 mic (SN  9945):   3.55 mV/Pa, -49.0 dB re 1V/Pa\n');
            case '7016#10252'
                fprintf(1, 'Calibration information from Manufacturer, 2018 (7016 #10252) \n');
                fprintf(1, '1/4\" 7016 mic (SN  10252):  3.85 mV/Pa, -48.3 dB re 1V/Pa\n');
        end
        fprintf(2, 'Transfer factor:         %7.2f mV/Pa, %5.1f dB re 1V/Pa\n', MIC.mVPerPa, MIC.dBPerVPa);
        fprintf(2, 'Measured with Mic Amp Gain = 20.0dB, Standard is %.1f dB\n', MIC.RefSig);
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
    fprintf(1, "Cosinor starts... ");
    for i = 1:length(fr)
        fprintf(1, 'computing i=%3d f=%7.3f\r', i, fr(i));
        [amp(i), theta(i), M(i), p_3a(i)] = cosinor(trec', ch2', 2.0*pi*fr(i), 0.05);
    end
    fprintf(1, "Cosinor done\n");

    return
end

function [amp_fft, fr_fft, x_pks, y_pks] = compute_fft(ch2, Fs, pad)
    Y = fft([detrend(ch2), zeros(1, length(ch2)*pad)]*(pad+1));  % only 1/5
    L = length(ch2)*(pad + 1);
    P2 = abs(Y/L);
    amp_fft = P2(1:L/2+1);
    amp_fft(2:end-1) = 2.0*amp_fft(2:end-1);
    fr_fft = Fs*(0:(L/2))/L;
    j = find(fr_fft >= 940 & fr_fft <= 1040);
    [pks, locs] = findpeaks(amp_fft(j));
    x_pks = fr_fft(j(1)-1+locs);
    y_pks = pks;
    return 
end
