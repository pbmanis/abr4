function [STIM] = updateStimParams(STIM, GUI)
% Get the stimulus parameters from the display window
% and store them in the STIM structure or vice versa,,
% probably need to clean this up... 


STIM.Info = 'ABR4 StimFile';
STIM.Date = date; % store date of last change

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # sweeps in average
set(GUI.hsr, 'string', num2str(STIM.NSweeps));
% Intersweep interval
set(GUI.hisi, 'String', num2str(STIM.InterSweepInterval));
% Stim/sweep
set(GUI.hsps, 'String', num2str(STIM.StimPerSweep));
% read interstimulus interval within a sweep
set(GUI.hipi, 'string', sprintf('%6.1f', STIM.ipi));
% Display duration
set(GUI.hdur, 'String', sprintf('%d', STIM.avg_dur));
% artificat rejection criterion: Amplitude criteria
set(GUI.hreject, 'String', num2str(STIM.avg_reject));
% artificat rejection criterion (stdev criteria)
STIM.rms_reject = str2double(get(GUI.hreject2, 'String'));
% amplifier gain setting - nominal is 10000.
set(GUI.hgain, 'String', sprintf('%d', STIM.amp_gain));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tone map duration. This is actually NOT saved...
%set(GUI.hToneMapDuration, 'string', STIM.tone_map_duration)
% lowest SPL in map or intensity series
set(GUI.hmin, 'String', num2str(STIM.SPLMin));
% highest SPL in map or intensity series
set(GUI.hmax, 'String', num2str(STIM.SPLMax));
% Intensity step size, dBSPL
set(GUI.hstep, 'String', num2str(STIM.SPLStep));
% Frequencies (kHz), as string list or shorhand
set(GUI.hflist, 'String', STIM.FreqList);
% set the alternation check box
set(GUI.hchk, 'value', STIM.Alternate);

% reset the rep counter
set(GUI.hrep, 'String', '0');
drawnow;