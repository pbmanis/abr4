function [varargout] = tonepip(amp, freq, delay, duration, rf, phase0, sampfreq, ipi, np, alternate)
% tonepip  - generate a tone pip with amplitude (V), frequency (Hz)
% delay (msec), duration (msec).
% if no rf (risefall) time is given, cosine shaping 5 msec is applied.
% if no phase is given, phase starts on 0, with positive slope.
%
% 4/2010 P. Manis. Minor fix for delay to first tone (was adding 2 delays).
% 

if nargout == 0 
   fprintf(2, 'amp: %f  freq: %f  delay = %f  duration = %f  rf = %f  ipi = %f  np = %f\n',...
       amp, freq, delay, duration, rf, ipi, np);
   return
end

if(nargin < 4)
    fprintf(2, 'Tonepip requires at least 4 arguments: amp, freq, delay, duration\n');
    return;
end
if(nargin == 4)
    phase0 = 0; % set defaults
    rf = 5;
    np = 1;
    ipi = 20;
    alternate  = 0;
end
if(nargin < 5)
    phase0 = 0;
    np = 1;
    ipi = 20;
    alternate  = 0;
end
if(nargin == 6)
    sampfreq = 500000; % sec....
    np = 1;
    ipi = 20;
    alternate  = 0;
end
% assume that ipi, np and alternate are set if more than 6 arguments.
if(nargin > 10)
    disp('too many arguments to tonepip');
    return;
end
if(nargout == 0)
    plotflag = 1;
else
    plotflag = 0;
end

clock = 1000/sampfreq; % calculate the sample clock rate - msec (khz)
uclock = 1000*clock; % microsecond clock
phi = 2*pi*phase0/360; % convert phase from degrees to radians...

jdur = floor((duration-2*rf)/clock); % duration of a signal
np_ipi = floor(ipi/clock); %
% build sin^2 rising and falling filter from 0 to 90 deg for shaping the waveform
%
nfilter_points = floor(rf/clock); % number of points in the filter rising/falling phase
fo = 1/(4*rf); % filter "frequency" in kHz - the 4 is because we use only 90deg for the rf component

fil = zeros((jdur+2*nfilter_points),1);

for i = 1:nfilter_points % rising filter shape
    fil(i) = sin(2*pi*fo*(i-1)*clock)^2; % filter
end
for j = i+1:i+jdur % main part shape
    fil(j) = 1;
end

i = 0;
for k = j+1:j+nfilter_points % decay shape
    fil(k) = fil(nfilter_points+i); %reverse the rising phase
    i = i - 1;
end
%Fs = 1000/clock;
%phi = 0; % initial phase
tfil = 0:clock:(length(fil)-1)*clock;
ws = amp*sin(phi + 2*pi*freq/1000*tfil)';
if rf > 0.0
    wf = ws.*fil; % this makes the stimulus pulse (sine, ramped)
else
    wf = ws;
end
nwf = length(wf);
%
% next put in context and make an output waveform
%
id = floor(ipi/clock); % spacing between pulses
jd = floor(delay/clock); % delay to start of stimulus

if(np > 1)
%    w = zeros(jd+np*np_ipi+nwf, 1);
    w = zeros(jd+np*np_ipi, 1);
else
    if jd + nwf < np_ipi
        w = zeros(np_ipi, 1);
    else
        w = zeros(jd+nwf, 1);
    end
end

for i = 1:np
    j0 = jd + (i-1)*id + 1;
    if jd == 0 && rf == 0.0
        j0 = 1 + (i-1)*id;
    end
    
    if(alternate && mod(i,2) == 1)
        sign = -1;
    else
        sign = -1;
    end
    ij=1;
    for j = j0:j0+nwf-1
        w(j) = sign*wf(ij);
        ij = ij + 1;
    end
    
end
w = w(1:length(w)-jd);  % cut tail points out of waveform so that we can
                        % concatenate without a delay.

if(nargout >= 1)
    varargout{1} = w;
end
if(nargout >= 2)
    varargout{2} = uclock;
end
if(nargout >= 3)
    varargout{3} = fil;
end

%fprintf(1, 'Stim Duration: %f pts,   %f msec\n', length(w), length(w)*clock);
if(plotflag)
    t = 0:clock:(length(w)-1)*clock;
    ff = findobj('tag', 'tonepip_figure');
    if isempty(ff)
        ff = figure;
        set(ff, 'tag', 'tonepip_figure');
        set(ff, 'Name', 'Tone Pip Test');
        set(ff, 'NumberTitle', 'off');
    else
        figure(ff);
        cla;
    end
    plot(t, w);
    hold on
    plot(tfil, fil, 'r');
    plot(tfil, ws, 'g');
    
    plot(tfil, wf, 'c');
end

return;