function [ output_args ] = SimpleAcq()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure(1);
clf;
RP=actxcontrol('rpco.x', [5 5 26 26]);
if(invoke(RP, 'connectrp2', 'usb', 1) == 0)
    error('failed to connect to rp2');
end;

% designate the sample frequency
% 0 = 6K, 1 = 12K, 2 = 25k, 3 = 50k, 4 = 100k, 5 = 200k, > 5 is not defined.
samp_cof_flag =4; % 4 is for 100 kHz
samp_flist = [6103.5256125, 122107.03125, 24414.0625, 48828.125, ...
    97656.25, 195312.5];
if(samp_cof_flag > 5)
    samp_cof_flag = 5;
end;
RP2COFFlag = samp_cof_flag;
sample_freq = samp_flist(samp_cof_flag+1);

TraceDuration = 1.0; % seconds

nRecordPoints = sample_freq/TraceDuration;

if RP.ClearCOF() == 0
    error('failed to clear cof');
end;
thisdir = pwd;
if (RP.LoadCOFsf([thisdir '\abrs.rcx'], RP2COFFlag) == 0)
    error ('failed to load abrs.rcx file');
end;
sfreq=RP.GetSFreq();
fprintf(1, 'true sample frequency: %.6f hz\n', sfreq);
fprintf(1, 'Points: %d\n', nRecordPoints);
sample_freq= sfreq;

%     fprintf(2, 'rp_setup: Failed to set tag value rec_size with n_record_points\n');
%     err = 1;
%     invoke(RP, 'halt');
% return;
% end;
RP.Run();
status = double(RP.GetStatus());
if bitget(double(status), 1) == 0;
    fprintf(2, 'rp_setup: Error connecting to RP2.1\n');
    err = 1;
   RP.halt;
return;
elseif bitget(double(status), 2) == 0;
    fprintf(2, 'rp_setup: Error loading circuit to RP2.1\n');
    err = 1;
   RP.halt;
return;
elseif bitget(double(status), 3) ==0
    fprintf(2, 'Error running circuit in RP2.1\n');
    err = 1;
   RP.halt;
return;
else
    % disp('circuit loaded and running');
end
RP.SetTagVal('REC_size', nRecordPoints*2);

if nargin > 1 && strcmp(varargin{1}, 'start')
    if RP.SoftTrg(1) == 0 % start.
        fprintf(2, 'failed to set trigger on RP2.1\n');
        err = 1;
   RP.halt;
return;
    end;
end
RP.SoftTrg(1); % start.

curindex=RP.GetTagVal('Index');
while(curindex  < nRecordPoints) % Checks to see if it has read into half the buffer
    curindex=RP.GetTagVal('Index');
    fprintf(2, 'curindex: %d\n', curindex);
end
RP.SoftTrg(2);
ch1=double(RP.ReadTagV( 'data_out1', 0, nRecordPoints));
ch2=double(RP.ReadTagV( 'data_out2', 0, nRecordPoints));

ch1(1:20)

RP.Halt;

plot(ch1, 'g-');
hold on;
plot(ch2, 'r-');

end



