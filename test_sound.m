function [ch2] = test_sound(freq, dbspl, duration)
% Generate a tone at the specified frequency and intensity. 
% if a 3rd argument is present, it will be the duration.
global STIM SPKR
if nargin <= 2
    duration = 5.0;
end

 notchfilter = designfilt('bandstopiir', 'FilterOrder', 20, ...
            'HalfPowerFrequency1', 39000.0, 'HalfPowerFrequency2', 41000.0, ...
            'SampleRate', STIM.sample_freq);
[Speaker, Mic] = getSpeakerMic();
CAL_L = load(sprintf('frequency_%s.cal', Speaker), '-mat'); % get calibration file
            CAL = CAL_L.CAL;
fprintf(1, 'Using Speaker Calibration from %s\n', CAL.Date);
MIC_L = load('microphone_7016#10252.cal', '-mat');
MIC = MIC_L.MIC;
disp(CAL.Freqs)
disp(CAL.dBSPL)
splatF = interp1(CAL.Freqs, CAL.dBSPL, freq, 'spline'); 
fprintf(1, "splatf: %7.2f", splatF);
[CAL.Freqs*1e-3; CAL.dBSPL']'
attn = splatF - dbspl + SPKR.attn;
disp(attn)
if(attn < 0)
    attn = 0.0;
end
nr = floor(STIM.sample_freq*duration);
[STIM.wave, STIM.clock] = tonepip(10.0, freq, 0., duration*1000., 1.0, 0, STIM.NIFreq, 10, 1, 0);
set_attn(attn);
[~, ch2, ~] = calstim(nr);
disp(rms(ch2));
ybp =  filter(notchfilter, ch2);
 disp(1e3*rms(ybp))
db = MIC.Gain + 20.*log10((1e3*rms(ybp)/3.85)/2e-5)
figure
plot(ybp);
% figure
% periodogram(ybp,rectwin(length(ybp)),length(ybp),STIM.sample_freq)

end

