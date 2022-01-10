function [w, sr] = noise_gen(varargin)
%
% noise generator for audtiory experiments
% This routine generates a noise using the input arguments as follows:
% varargin{1} is the mode
%     mode can be: wideband, bandpass, lowpass, highpass, octave
%
% The subsequent arguments depend on the mode
% wideband: no additional arguments
% bandpass: varargin{2} is one filter corner frequency, varargin{3} is the
% other. These are sorted later.
% lowpasss: varargin{2} is the corner frequencey
% highpass: varargin{2} is the corner frequency
% octave: varargin{2} is the center frequency, varargin{3} is the bandwith
% in octaves.
% The next two arguments, if given, are the filter type and the order.
% if not given, the filter is Butterworth, 8th order.
%
% Outputs:
% if None, then the program generates duration sec worth of data at 500 kHz, and
% calculates the power spectrum of the result
% otherwise,it fills the output with the waveform and the sample rate.

% 26 March 2007 Paul B. Manis, Ph.D.
% UNC Chapel Hill

% generate the wideband noise signal first
%
duration = 2.0; % second
Fs = 500000; % sample rate, points per second (Hz)
sr = 1.0/Fs;
npnts = floor(duration*Fs);
%npnts = 100000; % override...
wbn = 3*randn(1, npnts);
% splice start and end a bit (weighted)
% otherwise you hear a click
wbn = [wbn wbn];
%w0 = (wbn(1)+wbn(end))/2.0;
%wbn(end) = (wbn(1)+wbn(end))/2.0;
%wbn(1) = w0;
ftype = 'butter';
order = 6;
nstage = 1;

switch(varargin{1})
    case 'wideband'
        fclp = Fs;

        % do nothing
    case {'lowpass', 'LP', 'lp'}
        passtype = 'low';
        fclp = varargin{2};
        if(nargin > 2)
            ftype = varargin{3};
            if(nargin > 3)
                order = varargin{4};
            end;
            if(nargin > 4)
                nstage = varargin{5};
            end;
        end
    case {'highpass', 'HP', 'hp'}
        passtype = 'high';
        fclp = varargin{2};
        if(nargin > 2)
            ftype = varargin{3};
            if(nargin > 3)
                order = varargin{4};
            end
            if(nargin > 4)
                nstage = varargin{5};
            end;

        end
    case {'bandpass', 'BP', 'bp'}
        passtype = 'bandpass';
        [fclp] = sort([varargin{2} varargin{3}]);
        if(nargin > 3)
            ftype = varargin{4};
            if(nargin > 4)
                order = varargin{5};
            end;
            if(nargin > 5)
                nstage = varargin{6};
            end;
        end
    case {'notch', 'stop'}
        passtype = 'stop';
        [fclp] = sort([varargin{2} varargin{3}]);
        if(nargin > 3)
            ftype = varargin{4};
            if(nargin > 4)
                order = varargin{5};
            end;
            if(nargin > 5)
                nstage = varargin{5};
            end;
        end
    otherwise
end

w = wbn;
% cutoff frequency in Hz
wco = fclp/(Fs/2); % wco of 1 is for half of the sample rate, so set it like this...
if(all(wco < 1)) % if wco is > 1 then this is not a filter!
    [b, a] = filterselect(ftype, order, wco, passtype);
    for i = 1:nstage
        w = filter(b, a, w); % filter all the traces...repeatedly
    end;
    for i = 1:nstage
        [b, a] = filterselect('elliptic', order, wco, passtype);
        w = filter(b, a, w); % filter all the traces...repeatedly
    end
end
k = length(w)/4;
w = w(k:(length(w)/2+k)-1);
% hw = figure();
% plot([w(length(w)-50:end) w(1:50) ])
if(nargout == 0)
    hs = spectrum.yulear(1024);

    h = findobj('tag', 'noise_gen_fig');
    if(~isempty(h))
        figure(h);
        clf;
        psd(hs, w, 'Fs', Fs);
    else
        figure('tag', 'noise_gen_fig');
        psd(hs, w, 'Fs', Fs);
    end;
    set(gca, 'Xlim', [0.1 100]);
    set(gca, 'XScale', 'log');
end;

function [b, a] = filterselect(type, order, wco, passtype)

switch(type)
    case 'butter'
        [b, a] = butter(order, wco, passtype); % butterworth
    case 'elliptic'
        [b, a] = ellip(order, 0.1, 120, wco, passtype);
    case 'cheby1'
        [b, a] = cheby1(order, 0.1, wco, passtype);
    case 'cheby2'
        [b, a] = cheby2(order, 10, wco, passtype);
    otherwise
        [b, a] = butter(order, wco, passtype); % butterworth

end;
