function updateStimParams()
% Get the stimulus parameters from the display window
% and store them in the STIM structure
global STIM

STIM.Info = 'ABR4 StimFile';
STIM.Date = date; % store date of last change

% read the alternation check box
hchk = findobj('Tag', 'ABR_AlternatePolarity');
set(hchk, 'value', STIM.Alternate);

% # sweeps in average
hsr = findobj('Tag', 'ABR_NSweeps'); % number of sweeps to average across
set(hsr, 'string', num2str(STIM.NSweeps));
STIM.hsr = hsr; % need to keep the handle as well.

% rep counter
hrep = findobj('Tag', 'ABR_StimRep');
set(hrep, 'String', '0');
STIM.hrep = hrep;

hisi = findobj('Tag', 'ABR_InterSweepInterval'); % number of sweeps to average across
STIM.hisi = hisi; % need to keep the handle as well.
set(hisi, 'String', num2str(STIM.InterSweepInterval));

% read the number of STIMULI per sweep (high rate bursts)
hsps = findobj('Tag', 'ABR_StimPerSweep');
set(hsps, 'String', num2str(STIM.StimPerSweep));

% read interstimulus interval within a sweep
hipi = findobj('Tag', 'ABR_SweepIPI');
set(hipi, 'string', sprintf('%6.1f', STIM.ipi));

% Display duration
hdur = findobj('Tag', 'ABR_Duration');
set(hdur, 'String', sprintf('%d', STIM.avg_dur));

% artificat rejection criterion: AMplitude criteria
hreject = findobj('Tag', 'ABR_Reject');
set(hreject, 'String', num2str(STIM.avg_reject));

% artificat rejection criterion (stdev criteria)
hreject = findobj('Tag', 'ABR_Reject2');
STIM.rms_reject = str2double(get(hreject, 'String'));

% amplifier gain setting - nominal is 10000.
hgain = findobj('Tag', 'ABR_Gain');
set(hgain, 'String', sprintf('%d', STIM.amp_gain));

hmin = findobj('Tag', 'ABR_MinAttn');
set(hmin, 'String', num2str(STIM.SPLMin));

hmax = findobj('Tag', 'ABR_MaxAttn');
set(hmax, 'String', num2str(STIM.SPLMax));

hstep = findobj('Tag', 'ABR_AttnStep');
set(hstep, 'String', num2str(STIM.SPLStep));

fh = findobj('Tag', 'ABR_FreqList');
set(fh, 'String', STIM.FreqList);

drawnow;