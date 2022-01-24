function [STIM] = loadStim(STIM, GUI)

[FileName,PathName,~] = uigetfile('*.abr4','Stim File to Load', 'StimFiles/*.abr4');
if FileName == 0 % cancelled out.
    return;
end

s = load([PathName FileName], '-mat');
if isfield(s.STIM, 'Info') && strcmp(s.STIM.Info, 'ABR4 StimFile')
else
    fprintf(2, 'File does not appear to be an ABR4 Stim File\n');
    return;
end
% update our STIM structure with comparable fieldnames in the
% retrieved file (which may have a slightly different format if it
% is older than 2022).
fn_load = fieldnames(s.STIM);
fn_current = fieldnames(STIM);
for i = 1:numel(fn_load)
   % fprintf(1,'  %s\n', fn_old{i});
    tf = strcmp(fn_load{i}, fn_current);
    if sum(tf) == 1
        
        kf = find(tf, 1);
        STIM.(fn_current{kf}) = s.STIM.(fn_load{i})
%         disp('----------------------')
%         disp(kf)
%         % disp(fn_new{tf})
%         % disp(s.STIM.(fn_old{i}))
%          if ~ strcmp(fn_load{i}, 'wave')
%              disp(fn_load{i})
%          disp(s.STIM.(fn_load{i}))
%          disp('===')
%          disp(fn_current{kf})
%          STIM.(fn_current{kf})
%          disp(STIM.(fn_current{kf}))
%          disp('***************************')
%          end  
    end
end
STIM = updateStimParams(STIM, GUI); % store info back to the window...
if ~isempty(GUI.hstimfilename)
    set(GUI.hstimfilename, 'String', FileName);
end
end
