function [h] = get_axis(PLOT, varargin)

arg = varargin{1};
if(nargin > 1)
    erflag = varargin{2}; % pick up erase flag
else
    erflag = 0; % no erase...
end
switch(arg)
    case {'signal', 'signal1'} 
        h = PLOT.signal1;
    case 'signal2'
        h = PLOT.signal2;
    case 'data'
        h = PLOT.Pdata;
    case 'response'
        h = PLOT.responsemap;
    otherwise
        h = [];
        error('? bad call to get_axis - arg = %s\n', arg);
end

