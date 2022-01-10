function [h] = get_axis(varargin)
global PLOTHANDLES

arg = varargin{1};
if(nargin > 1)
    erflag = varargin{2}; % pick up erase flag
else
    erflag = 0; % no erase...
end;
switch(arg)
    case {'signal', 'signal1'} 
        h = PLOTHANDLES.signal1;
    case 'signal2'
        h = PLOTHANDLES.signal2;
    case 'data'
        h = PLOTHANDLES.data;
    case 'response'
        h = PLOTHANDLES.responsemap;
    otherwise
        h = [];
        error('? bad call to get_axis - arg = %s\n', arg);
        return;
end;
 return
 
if(erflag && ~isempty(h))
    cla(h);
    hold off % force clear on next draw
end;

