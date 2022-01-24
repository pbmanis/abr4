function clear_plots(varargin)
% just clear the data and signal plots
global PLOTHANDLES REFERENCE ACQPars
cla(PLOTHANDLES.data);
cla(PLOTHANDLES.signal1);
cla(PLOTHANDLES.signal2);
set(PLOTHANDLES.data, 'YLimMode', 'auto');
set_scale();
if nargin > 0
    if ~isempty(REFERENCE)
        plot(PLOTHANDLES.data,  ACQPars.tb, REFERENCE, 'k-', 'linewidth', 2.0);
    end
end

return;