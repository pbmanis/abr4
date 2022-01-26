function click_cal()
%Click calibration:  Calibrate a click against a tone,
% using peak-to-peak method (Burkard, 2006).


H = abr4('gethw');
[Speaker, Mic] = getSpeakerMic();
MICL = load(sprintf('microphone_%s.cal', Mic), '-mat'); % get the current microphone calibration file
MIC = MICL.MIC;
STIM = H.STIM;
SPKR = H.CALIBRATION.SPKR;
HW = H.HW;
click_attn = 10.0;

[ch2, p2p] = test_sound(16001, 75, 1, false);
recordDuration = 1.0;
STIM.rate = 1.0/STIM.sample_freq;
STIM.click_dur = 0.1;  % msec
STIM.click_delay = 0.5;  % msec
STIM.NSweeps = 1;
STIM.ipi = 10;  % msec
STIM.Alternate = 0;
STIM.StimPerSweep = floor(1000*recordDuration/STIM.ipi)-1;

if strcmp(HW.HARDWARE, 'NI')
    STIM.NIFreq = 500000; % express in sec per point
else
    STIM.NIFreq = 44100;
end
STIM.delay = STIM.click_delay;
%STIM = updateStimParams(STIM, GUI);
fprintf(2, "Clicks: ");
HW = set_attn(HW, click_attn);
% click(amp, delay, duration, samplefreq, ipi, np, alternate)
[STIM.wave,STIM.clock] = click(H.CALIBRATION.SPLCAL.click_amp, STIM.click_delay, STIM.click_dur,...
    STIM.NIFreq, STIM.ipi, STIM.StimPerSweep, STIM.Alternate); % convert rate to usec per point
nRecordPoints = floor(recordDuration*STIM.sample_freq);
length(STIM.wave);
[~, ch2, HW, ~] = calstim(nRecordPoints, HW, STIM);
% design a notch filter at 40.0 kHz
notchfilter1 = designfilt('bandstopiir', 'FilterOrder', 20, ...
    'HalfPowerFrequency1', 39000.0, 'HalfPowerFrequency2', 41000.0, ...
    'SampleRate', STIM.sample_freq);
ch2 = filter(notchfilter1, ch2);
for ncenter = [800, 2400, 5200]
    notchfilter2 = designfilt('bandstopiir', 'FilterOrder', 20, ...
        'HalfPowerFrequency1', ncenter-50, 'HalfPowerFrequency2', ncenter+50, ...
        'SampleRate', STIM.sample_freq);
    ch2 = filter(notchfilter2, ch2);
end
tb = (0: 1.0/STIM.sample_freq: (length(ch2)-1)/STIM.sample_freq);
sttb = (0: 1.0/STIM.clock: (length(STIM.wave)-1)/STIM.clock);
HW = set_attn(HW, 0);
figure(97)
subplot(3,1,1);
plot(tb, ch2);
subplot(3,1,2);
plot(sttb, STIM.wave);

ts1 = floor(0.3*STIM.sample_freq);
ts2 = floor(0.95*STIM.sample_freq);
ch_c = ch2(ts1:ts2);
i_rate = 2e-6;
sps = floor((0.95-0.3)/(1e-3*STIM.ipi));
tb_interpolate = (tb(ts1): i_rate : tb(ts2));
min(tb_interpolate);
max(tb_interpolate);
ch_i = interp1(tb(ts1:ts2), ch_c,tb_interpolate, 'linear', 'extrap');

interBlockLen = floor(1e-3*STIM.ipi/i_rate);
fprintf(2, "interblock length: %d   vs. %8.3f", interBlockLen, 1e-3*STIM.ipi/i_rate)
% length(ch_i)
% interBlockLen
% sps
% interBlockLen*sps
if interBlockLen*sps > length(ch_i)
    sps = floor(length(ch_i)/interBlockLen);
end

ch_r = reshape(ch_i(1:interBlockLen*sps), ...
    [interBlockLen, sps]);
% size(ch_r)
ch_rm = mean(ch_r, 2);
subplot(3,1,3);

tb_r = (0:i_rate:(length(ch_rm)-1)*i_rate);
%plot(tb_r, ch_rm);
%[bp_data, bp_freqs] = periodogram(ch_rm,rectwin(length(ch_rm)),length(ch_rm),1/i_rate);
%plot(bp_freqs, bp_data);
size(ch_r)
for i = 1:size(ch_r, 2)
    if i % 5 ~= 0
        continue
    end
    plot(tb_r*1e3, ch_r(:, i));
    hold on
end
plot(tb_r, ch_rm, 'k-');
p2p = peak2peak(ch_rm);
clickspl = soundfuncs.compute_spl(p2p/(2*sqrt(2)), MIC);  % based on RMS
fprintf(2, "Peak to Peak click: %8.3f mV\n",p2p*1e3);
fprintf(2, "dBSPL = %8.2f dB with %5.1f dB Attn, maxclick = %5.1f dBSPL\n", ...
    clickspl, click_attn, clickspl+click_attn);
end

