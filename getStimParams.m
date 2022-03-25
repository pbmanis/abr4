function [STIM] = getStimParams(STIM, GUI)
% Get the stimulus parameters from the display window
% and store them in the STIM structure

STIM.Info = 'ABR4 StimFile';
STIM.period = 100;
STIM.risefall = 0.1;
[~, STIM] = get_spls(STIM, GUI);
[~, STIM] = get_freqs(STIM, GUI);
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
STIM.NSweeps = str2double(get(GUI.hsr, 'String'));

% display rep counter
set(GUI.hrep, 'String', '0');
STIM.InterSweepInterval = str2double(get(GUI.hisi, 'String'));

% read the number of STIMULI per sweep (high rate bursts)
hsps = findobj('Tag', 'ABR_StimPerSweep');
if(isempty(GUI.hsps))
    STIM.StimPerSweep = 20;
else
    sps = str2double(get(GUI.hsps, 'String'));
    if(sps < 1 || sps > 100)
        sps = 20;
        set(hsps, 'string', sprintf('%d', sps));
    end
    STIM.StimPerSweep = sps;
end

% read interstimulus interval within a sweep
hipi = findobj('Tag', 'ABR_SweepIPI');
if(isempty(GUI.hipi))
    STIM.ipi = 40;
else
    ipi = str2double(get(GUI.hipi, 'String'));
    if(ipi < 1 || ipi > 5000)
        ipi = 40;
        set(GUI.hipi, 'string', sprintf('%6.1f', ipi));
    end
    STIM.ipi = ipi;
end

% Display duration
STIM.avg_dur = str2double(get(GUI.hdur, 'String'));
if(STIM.avg_dur <= 5 || STIM.avg_dur > STIM.ipi)
    STIM.avg_dur = STIM.ipi; % milliseconds
    set(GUI.hdur, 'String', sprintf('%d', STIM.avg_dur));
end

% artificat rejection criterion: Amplitude criteria
STIM.avg_reject = str2double(get(GUI.hreject, 'String'));

% artificat rejection criterion (stdev criteria)
STIM.rms_reject = str2double(get(GUI.hreject2, 'String'));

% amplifier gain setting - nominal is 10000.
amp_gain = str2double(get(GUI.hgain, 'String'));
if(amp_gain < 1 || amp_gain > 1000000)
    amp_gain = 10000; % gain
    set(GUI.hgain, 'String', sprintf('%d', amp_gain));
end
STIM.amp_gain = amp_gain;


function [spls, STIM] = get_spls(STIM, GUI)
% return the sound pressure level sequence from the window
%
% the tags in the windows are attenuations, but in this version
% we treat as if SPL
%
minspl = str2double(get(GUI.hmin, 'string'));
maxspl = str2double(get(GUI.hmax, 'string'));
% hstep = findobj('Tag', 'ABR_AttnStep');
stepspl = str2double(get(GUI.hstep, 'string'));
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

function [freqs, STIM] = get_freqs(STIM, GUI)
% return the frequency list from the window
fllist = get(GUI.hflist, 'string');
STIM.FreqList = fllist; % save the actual list as well.
freqs = seq_parse(fllist);
if(max([freqs{:}]) < 100) % then they meant kHz, not Hz !
    for i = 1:length(freqs)
        freqs{i} = freqs{i}*1000;
    end
end
STIM.freqs = cell2mat(freqs);
return;
