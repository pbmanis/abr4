function datac_cursorread(n)
%
% get the most recent cursor click position 
% also redraws/creates the updated marker lines if n < 0...
global TMARKS % access to the time markers on the display.

h_axes = findobj('tag', 'datawin1');
if(isempty(h_axes) | ~ishandle(h_axes))
   return;
end;
axinfo = get(h_axes, 'UserData'); % There may be information in userdata we can use
if(isempty(axinfo) | isfield(axinfo, 'PlotHandle')) % not OUR baby...
   return;
end;
if(ischar(n))
   n = str2num(n);
end;
if(n >= 0)
%   disp 'updatemarks'
   htgt = findobj('tag', sprintf('data_text%1d', n));
   if(~isempty(htgt) & ishandle(htgt))
      set(htgt, 'String', sprintf('%7.2f', axinfo.down_pos(1)));
      % fprintf(1,'%7.2f\n', axinfo.down_pos(1));
      set(htgt, 'value', axinfo.down_pos);
      TMARKS.t(n+1) = axinfo.down_pos(1);
      if(n+1 < length(TMARKS.t))
         set(TMARKS.h(n+1), 'XData', [TMARKS.t(n+1) TMARKS.t(n+1)]);
      	set(TMARKS.h(n+1), 'erasemode', 'xor');
   end;
end;
else % we draw the cursor lines on the graph
   %disp 'drawcursorlines'
   set(gcf, 'CurrentAxes', findobj('tag', 'datawin1'))
   if(length(TMARKS.t) < 4)
      ntm = length(TMARKS.t);
   else
      ntm = 4;
   end;
   for i = 1:ntm
      hold on;
      u = get(gca, 'Ylim');
      if(~isfield(TMARKS, 'h')) % no line handles - need to create the lines
         TMARKS.h(i)=plot([TMARKS.t(i) TMARKS.t(i)], u, 'r-', 'erasemode', 'xor');
      else % field exists - see if/how we can use it
         if(i <= length(TMARKS.h)) % array exists - and is long enough 
            if(~isempty(TMARKS.h(i)) & ishandle(TMARKS.h(i))) % fill it 
               set(TMARKS.h(i), 'XData', [TMARKS.t(i) TMARKS.t(i)]);
            else % not a handle or maybe empty - make one...
               TMARKS.h(i)=plot([TMARKS.t(i) TMARKS.t(i)], u, 'r-', 'erasemode', 'xor');
				end;
         else % not enough existing - make a new one
            TMARKS.h(i)=plot([TMARKS.t(i) TMARKS.t(i)], u, 'r-', 'erasemode', 'xor');
         end;      
      end;      
      
   end;
   
end;

