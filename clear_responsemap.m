function clear_responsemap() 
global RESPONSE_MAP
get_axis('response');
if(isempty(RESPONSE_MAP) || ~ishandle(RESPONSE_MAP))
    RESPONSE_MAP = [];
    ylim('auto');
    RESPONSE_MAP(1) = plot([0 0], [0 0]);
    set(RESPONSE_MAP(1), 'linestyle', '-', 'marker', 'x', 'color', 'green');
else
    ylim('auto');
    set(RESPONSE_MAP(1), 'Xdata', [0 0], 'Ydata', [0 0]);
end;
drawnow;
end