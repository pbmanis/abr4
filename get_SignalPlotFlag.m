function [flag] = get_SignalPlotFlag(signo)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 h = findobj('Tag', sprintf('ABR_Stim%d_Show', signo));
 if ~isempty(h)
     flag = get(h, 'Value');
 else
     flag = 0;
 end;
 
end

