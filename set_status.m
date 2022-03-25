function set_status( statustext)
%%
% set status indicator
h = findobj('Tag', 'ABR_Status');
set(h, 'String', statustext);
h = findobj('Tag', 'ABR_Stopbutton');
set(h, 'UserData', statustext);
end

