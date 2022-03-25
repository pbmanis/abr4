function [state] = check_status(GUI)
% get the program state from the GUI stop button
    state = get(GUI.hstop, 'UserData');
end


