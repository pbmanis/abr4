function clear_responsemap(responsemap) 

get_axis('response');
if(isempty(responsemap) || ~ishandle(responsemap))
    responsemap = [];
    ylim('auto');
    responsemap(1) = plot([0 0], [0 0]);
    set(responsemap(1), 'linestyle', '-', 'marker', 'x', 'color', 'green');
else
    ylim('auto');
    set(responsemap(1), 'Xdata', [0 0], 'Ydata', [0 0]);
end
drawnow;
end