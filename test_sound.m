function [ch2, p2p] = test_sound(freq, dbspl, duration, display)
% Generate a tone at the specified frequency and intensity.
% if a 3rd argument is present, it will be the duration.

H = abr4('gethw');

STIM = H.STIM;
SPKR = H.CALIBRATION.SPKR;
HW = H.HW;

if nargin <= 2
    duration = 5.0;
end
if nargin <= 3
    display = true;
end

notchfilter = designfilt('bandstopiir', 'FilterOrder', 20, ...
    'HalfPowerFrequency1', 39000.0, 'HalfPowerFrequency2', 41000.0, ...
    'SampleRate', STIM.sample_freq);
[Speaker, Mic] = getSpeakerMic();

CAL_L = load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file
CAL = CAL_L.CAL;
fprintf(1, 'Using Speaker Calibration from %s\n', CAL.Date);
micfile = sprintf('microphone_%s.cal', Mic);
fprintf(1, 'Using microphone calibration from:%s\n', micfile);
MIC_L = load(micfile, '-mat');
MIC = MIC_L.MIC;
% disp(CAL);
% disp(MIC);
splatF = soundfuncs.spl_at_f(CAL.Freqs, CAL.maxdB, freq);
% interp1(CAL.Freqs, CAL.dBSPL, freq, 'spline');
fprintf(1, "splatf: %7.2f, freq: %7.1f\n", splatF, freq);
[CAL.Freqs*1e-3; CAL.dBSPL']';
attn = splatF - dbspl;
% figure
% plot(CAL.Freqs, CAL.dBSPL, 'ko');
% spl_freqs = (min(CAL.Freqs):100:max(CAL.Freqs));
% splatF2 = interp1(CAL.Freqs, CAL.dBSPL, spl_freqs, 'spline');
% hold on
% plot(spl_freqs, splatF2, 'r-');
if(attn < 0)
    attn = 0.0;
end
nr = floor(STIM.sample_freq*duration*1.2);
[STIM.wave, STIM.clock] = tonepip(10.0, freq, 0, duration*1000., 10.0, 0, STIM.NIFreq, 10, 1, 0);
set_attn(HW, attn);
fprintf(1, "Attenuator: %5.1f\n", attn);
tstim = (0:1/STIM.NIFreq:(length(STIM.wave)-1)/STIM.NIFreq);
[~, ch2, ~] = calstim(nr, HW, STIM);
ybp =  ch2; % filter(notchfilter, ch2);
trec = (0:1/STIM.sample_freq:(nr-1)/STIM.sample_freq);
ts1 = floor(0.1*duration*STIM.sample_freq);
ts2 = floor(0.9*duration*STIM.sample_freq);
fprintf(1, "Raw Vrms:       %8.3f mV\n",(1e3*rms(ybp)));
[amp_cosinor, fr_cosinor, ~, ~] = cosinor(trec(ts1:ts2)', ybp(ts1:ts2)', 2.0*pi*freq, 0.05);

amp_cosinor_rms_V = amp_cosinor/sqrt(2);
bpfreqs = soundfuncs.octave_calc(freq, 5, STIM.sample_freq);
%             fprintf(1, "spkr: %f  bp: %f  to %f\n", spkr_freq(i), bpfreqs(1), bpfreqs(2));
ybp = bandpass(ch2, bpfreqs, STIM.sample_freq, ...
    'StopbandAttenuation', 60, "Steepness", 0.9);
p2p = sqrt(2)*2*rms(ybp);
fprintf(1, "Signal Vrms with Gain = %.1f dB Vsig_cosinor = %8.2f mV, bandpass: %8.2f mV Mic Vref: %8.5f V\n", ...
    MIC.Gain, amp_cosinor_rms_V*1e3, rms(ybp)*1e3, MIC.Vref_bp);
db = soundfuncs.compute_spl(amp_cosinor_rms_V, MIC);
fprintf(1, "Measured microphone voltage P2P  %7.3f mV at %7.1f dBSPL:\n", p2p*1e3, db);
if display
    figure(98);
    subplot(2,1,1);
    plot(trec(1:ts1), ybp(1:ts1), 'r-');
    hold on
    plot(trec(ts2:end), ybp(ts2:end), 'r-');
    plot(trec(ts1:ts2), ybp(ts1:ts2), 'k-');
    subplot(2,1,2);
    plot(tstim, STIM.wave);
end

end

