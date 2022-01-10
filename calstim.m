function [ch1, ch2, err] = calstim(nRecordPoints)
%% Load and record signal to do a calibration 
%% assumes stimulus waveform is already in STIM.wave and that 
%% AO and RP are already setup.
%% 5/1/2010 Paul B. Manis
%% converted to matlab session interface to NIDAQ 9/28/2016
%%
global AO STIM RP

% fprintf(1, 'Stim clock: %8.1f usec, NIFreq: %8.1f Hz, Pts: %d\n', ...
%     STIM.clock, STIM.NIFreq, length(STIM.wave));

AO.stop 
AO.Rate = STIM.NIFreq; 
% queueOutputData(AO, STIM.wave); % wave is FULL
AO.TriggersPerRun = 1;
queueOutputData(AO, STIM.wave); % wave is FULL

if isempty(RP)
    RP=actxcontrol('rpco.x', [5 5 26 26]);
    if(RP.ConnectRP2('USB',1) == 0)
        error('Failed to connect to RP2.1');
    end;
end

% disp 'rp setup'
err = rp_setup(nRecordPoints, 'microphone');
startBackground(AO); % get ni board read to go, then trigger the rp
% pause(0.001);
RP.SoftTrg(1); % start.
% disp 'started'

%nRecordPoints = floor(STIM.sample_freq*TraceDuration*1.5)
TraceDuration = nRecordPoints/STIM.sample_freq;
curindex=RP.GetTagVal('Index');
lastindex = curindex;
% while(curindex  < nRecordPoints) % Checks to see if it has read into half the buffer
%     pause(0.025);
%     curindex=RP.GetTagVal('Index');
%     if curindex < lastindex
%         break
%     end;
%     lastindex = curindex;
% 
% 
% end
pause(TraceDuration)
RP.SoftTrg(2);
% disp 'stopped'
% ch1=double(RP.ReadTagV('data_out1', 0, nRecordPoints));
ch1 = 0;
ch2=double(RP.ReadTagV('data_out2', 0, nRecordPoints));
AO.stop;
set_attn(-1);
RP.SoftTrg(2);
RP.Halt;
end       