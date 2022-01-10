function datac_setcrosshair(axis_handle, tag, xu, yu, pos, callback, delta, type)
%
% datac_setcrosshair; given an axis and a named tag for the text window to
% display the measurements, set up for a crosshair callback for the axis.
%

% check the arguments.

if(nargin < 5)
   error('data_setcrosshair has insufficient arguments  in call');
end;
if(~ishandle(axis_handle))
   error('data_setcrosshair: axis_handle argument is not correct');
end;
if(~ischar(tag))
   error('data_setcrosshair: tag argument is bad');
end;
if(~ischar(xu))
   error('data_setcrosshair: xunit argument is bad');
end;
if(~ischar(yu))
   error('data_setcrosshair: yunit argument is bad');
end;
if(~isnumeric(pos) & length(pos) ~= 4)
   error('datac_setcrosshair: position for display field is bad');
end;
if(nargin > 6)
   if(~isnumeric(delta))
      error('datac_setcrosshair: delta flag must be numeric');
   end;
else
   delta = 0;
end;
if(nargin < 7)
   type = 'rectangle';
end;
switch(type)
case {'rectangle', 'time_marker'}
otherwise
   error('datac_setcrosshair: crosshair type not recognized');
   return;
end;

% safe to proceed:

htext = findobj('Tag', tag); % find/create the display window information
if(isempty(htext))
   if(delta)
      htext = uicontrol(...
      'units', 'normalized', ...
      'Position', pos, ...
      'Style', 'text', ...
      'string', sprintf(' %7.2f  %s\n %7.2f  %s\n*%7.2f   %s\n*%7.2f    %s', 0, xu, 0, yu, 0, xu, 0, yu),  ...
      'Tag', tag );
   else
      htext = uicontrol(...
      'units', 'normalized', ...
      'Position', pos, ...
      'Style', 'text', ...
      'string', sprintf('%7.1f  %s\n%7.1f  %s', 0, xu, 0, yu),  ...
      'Tag', tag );
   end;
else
   set(htext, 'Visible', 'on');
end;

% build the axinfo structure.

axinfo.text_handle = htext;
axinfo.type = type;
axinfo.text_position = pos;
axinfo.delta_mode = 0;
axinfo.delta = [0, 0];
axinfo.line_handle = [];
axinfo.xunit = xu;
axinfo.yunit = yu;
axinfo.curpos = [0 0]; % current crosshair position (absolute in axis coordinates)
if(nargin > 5)
   axinfo.callback = callback;
end;

% store the structure.
if(~isempty(axis_handle))
   set(axis_handle, 'UserData', axinfo);
end;

return;
