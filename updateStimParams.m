function [STIM] = updateStimParams(STIM, GUI)
% Get the stimulus parameters from the display window
% and store them in the STIM structure or vice versa,,
% probably need to clean this up... 

STIM.Info = 'ABR4 StimFile';
STIM.Date = date; % store date of last change

% read the alternation check box
set(GUI.hchk, 'value', STIM.Alternate);

% # sweeps in average
set(GUI.hsr, 'string', num2str(STIM.NSweeps));
STIM.hsr = hsr; % need to keep the handle as well.

% rep counter
set(GUI.hrep, 'String', '0');
STIM.hrep = hrep;

STIM.hisi = hisi; % need to keep the handle as well.
set(GUI.hisi, 'String', num2str(STIM.InterSweepInterval));

% read the number of STIMULI per sweep (high rate bursts)
set(GUI.hsps, 'String', num2str(STIM.StimPerSweep));

% read interstimulus interval within a sweep
set(GUI.hipi, 'string', sprintf('%6.1f', STIM.ipi));

% Display duration
set(GUI.hdur, 'String', sprintf('%d', STIM.avg_dur));

% artificat rejection criterion: AMplitude criteria
set(GUI.hreject, 'String', num2str(STIM.avg_reject));

% artificat rejection criterion (stdev criteria)
STIM.rms_reject = str2double(get(GUI.hreject2, 'String'));

% amplifier gain setting - nominal is 10000.
set(GUI.hgain, 'String', sprintf('%d', STIM.amp_gain));

set(GUI.hmin, 'String', num2str(STIM.SPLMin));

set(GUI.hmax, 'String', num2str(STIM.SPLMax));

set(GUI.hstep, 'String', num2str(STIM.SPLStep));

set(GUI.hflist, 'String', STIM.FreqList);

drawnow;