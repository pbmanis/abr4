function [varargout] = click(amp, delay, duration, samplefreq, ipi, np, alternate)
%
% compute click waveform or train with specified amplitude delay and duration
% and samplefreq (in Hz) (required)
% or train of identical clicks with ipi and np
% all time arguments in msec; samplefreq is in Hz.
% returns w (waveform) and uclock (max dac sample rate, in microseconds).
%
% 4/2010 : minor fix for the delay (was not calculated correctly). 

w=[];
clock = [];
if(nargin < 4)
    disp('click requires at least 3 arguments: amp, delay, duration and rate');
    return;
end;
if(nargin ~= 4 && nargin ~= 6 && nargin ~= 7)
    disp('click requires 4 (single click) or 6 (click train with ipi and np) arguments');
    return;
end;
if(nargin ==4 )
    np = 1;
    ipi = 1;
end;
if(nargin < 6)
    alternate = 0;
end;
if(nargout == 0)
    plotflag = 1;
else
    plotflag = 0;
end;

clock = 1000/samplefreq; % convert Hz to milliseconds.
uclock = 1000/clock; % convert to microseconds.
% gcd(1000*delay, 1000*duration)/1000;
%maxt = delay+ipi*np+2;
maxt = ipi*np;
% maxt = 10; % in milliseconds?
w = zeros(floor(maxt/clock), 1);
if(nargout == 0)
    fprintf(1, 'clock: %f   maxt: %f  array: %d\n', clock, maxt, length(w));
end;

jd = floor(delay/clock)+1; % delay time to first pulse in clock cycles
pd = floor(duration/clock); % pulse duration in clock cycles
id = floor(ipi/clock); % interpulse interval in clock cycles
for i = 1:np
    j0 = (i-1)*id + jd;
    j1 = j0 + pd;
    if(alternate && mod(i,2) == 1)
        sign = -1;
    else
        sign = 1;
    end;
    for j = j0:j1
        w(j) = sign*amp;
    end;
    
end;
if(nargout > 0)
    varargout{1} = w;
    varargout{2} = uclock;
end;
if(plotflag)
    t = [0:clock:(length(w)-1)*clock];
    figure;
    plot(t, w);
end;
