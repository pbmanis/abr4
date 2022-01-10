function set_status( statustext)
%%
% set status indicator
h = findobj('Tag', 'ABR_Status');
set(h, 'String', statustext);

end

