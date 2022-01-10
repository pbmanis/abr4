function set_scale()
%read the display scale list and set the plot scale
%
global PLOTHANDLES
ds = findobj('Tag', 'ABR_DisplayScale');
dsstr = get(ds, 'String');
dsval = get(ds, 'Value');
if strcmp(dsstr{dsval}, 'auto') == 0
    dsscale = str2double(dsstr{dsval})*1e-6;
    set(PLOTHANDLES.data, 'YLim', [-dsscale dsscale]);
else
        set(PLOTHANDLES.data, 'yLimMode', 'Auto');
end;
end

