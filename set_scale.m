function set_scale(PLOT)
%read the display scale list and set the plot scale
%

ds = findobj('Tag', 'ABR_DisplayScale');
dsstr = get(ds, 'String');
dsval = get(ds, 'Value');
if strcmp(dsstr{dsval}, 'auto') == 0
    dsscale = str2double(dsstr{dsval})*1e-6;
    set(PLOT.data, 'YLim', [-dsscale dsscale]);
else
        set(PLOT.data, 'yLimMode', 'Auto');
end
end

