function [STIM] = getStimParams(STIM)
% Get the stimulus parameters from the display window
% and store them in the STIM structure


STIM.Info = 'ABR4 StimFile';
STIM.period = 100;
STIM.risefall = 0.1;
[spls, STIM] = get_spls(STIM);
[freqs, STIM] = get_freqs(STIM);
STIM.NIFreq = 500000;
STIM.tone_delay = 1.0;
STIM.tone_delay_mapping = 10.0;
STIM.click_delay = 1.0;
STIM.StimPerSweep  = 1;

% read the alternation check box
hchk = findobj('Tag', 'ABR_AlternatePolarity');
if(ishandle(hchk))
    c = get(hchk, 'value');
    
    if(c == 1)
        STIM.Alternate = 1;
    else
        STIM.Alternate = 0;
    end
else
    STIM.Alternate = 0;
end

% get # sweeps in average
hsr = findobj('Tag', 'ABR_NSweeps'); % number of sweeps to average across
STIM.hsr = hsr; % need to keep the handle as well.
STIM.NSweeps = str2double(get(hsr, 'String'));

% display rep counter
hrep = findobj('Tag', 'ABR_StimRep');
set(hrep, 'String', '0');
STIM.hrep = hrep;

hisi = findobj('Tag', 'ABR_InterSweepInterval'); % number of sweeps to average across
STIM.hisi = hisi; % need to keep the handle as well.
STIM.InterSweepInterval = str2double(get(hisi, 'String'));

% read the number of STIMULI per sweep (high rate bursts)
hsps = findobj('Tag', 'ABR_StimPerSweep');
if(isempty(hsps))
    STIM.StimPerSweep = 20;
else
    sps = str2double(get(hsps, 'String'));
    if(sps < 1 || sps > 100)
        sps = 20;
        set(hsps, 'string', sprintf('%d', sps));
    end
    STIM.StimPerSweep = sps;
end

% read interstimulus interval within a sweep
hipi = findobj('Tag', 'ABR_SweepIPI');
if(isempty(hipi))
    STIM.ipi = 40;
else
    ipi = str2double(get(hipi, 'String'));
    if(ipi < 1 || ipi > 5000)
        ipi = 40;
        set(hipi, 'string', sprintf('%6.1f', ipi));
    end
    STIM.ipi = ipi;
end

% Display duration
hdur = findobj('Tag', 'ABR_Duration');
STIM.avg_dur = str2double(get(hdur, 'String'));
if(STIM.avg_dur <= 5 || STIM.avg_dur > STIM.ipi)
    STIM.avg_dur = STIM.ipi; % milliseconds
    set(hdur, 'String', sprintf('%d', STIM.avg_dur));
end

% artificat rejection criterion: Amplitude criteria
hreject = findobj('Tag', 'ABR_Reject');
STIM.avg_reject = str2double(get(hreject, 'String'));

% artificat rejection criterion (stdev criteria)
hreject = findobj('Tag', 'ABR_Reject2');
STIM.rms_reject = str2double(get(hreject, 'String'));

% amplifier gain setting - nominal is 10000.
hgain = findobj('Tag', 'ABR_Gain');
amp_gain = str2double(get(hgain, 'String'));
if(amp_gain < 1 || amp_gain > 1000000)
    amp_gain = 10000; % gain
    set(hgain, 'String', sprintf('%d', amp_gain));
end
STIM.amp_gain = amp_gain;


function [spls, STIM] = get_spls(STIM)
% return the sound pressure level sequence from the window
%

minspl = 0;
maxspl = 120;
stepspl = 10;
% the tags in the windows are attenuations, but in this version
% we treat as if SPL
%
hmin = findobj('Tag', 'ABR_MinAttn');
if(~isempty(hmin))
    minspl = str2double(get(hmin, 'string'));
end
hmax = findobj('Tag', 'ABR_MaxAttn');
if(~isempty(hmax))
    maxspl = str2double(get(hmax, 'string'));
end
hstep = findobj('Tag', 'ABR_AttnStep');
if(~isempty(hstep))
    stepspl = str2double(get(hstep, 'string'));
end
if(minspl > maxspl)
    temp = maxspl;
    maxspl = minspl;
    minspl = temp;
end
spls = minspl:stepspl:maxspl;
STIM.spls = spls;
STIM.SPLMin = minspl;
STIM.SPLMax = maxspl;
STIM.SPLStep = stepspl;
return;

function [freqs, STIM] = get_freqs(STIM)
% return the frequency list from the window

fh = findobj('Tag', 'ABR_FreqList');
fllist = get(fh, 'string');
STIM.FreqList = fllist; % save the actual list as well.
freqs = seq_parse(fllist);
if(max([freqs{:}]) < 100) % then they meant kHz, not Hz !
    for i = 1:length(freqs)
        freqs{i} = freqs{i}*1000;
    end
end
STIM.FreqList = freqs;
return;
