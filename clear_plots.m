function clear_plots(PLOTS, STIM, DATA, varargin)
% just clear the data and signal plots
cla(PLOTS.data);
cla(PLOTS.signal1);
cla(PLOTS.signal2);
set(PLOTS.data, 'YLimMode', 'auto');
set_scale(PLOTS);

if nargin > 1
    if ~isempty(DATA.REFERENCE) && (length(STIM.ACQPars_tb) == length(DATA.REFERENCE))
        plot(PLOTS.data,  STIM.ACQPars_tb, DATA.REFERENCE, 'k-', 'linewidth', 2.0);
    end
end

return;