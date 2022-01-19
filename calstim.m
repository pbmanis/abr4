function [ch1, ch2, err] = calstim(duration)
% Load and record signal to do a calibration 
% assumes stimulus waveform is already in STIM.wave and that 
% AO and RP are already setup.
% 5/1/2010 Paul B. Manis
% converted to matlab session interface to NIDAQ 9/28/2016
%
global AO RP STIM

% fprintf(1, 'Stim clock: %8.1f usec, NIFreq: %8.1f Hz, Pts: %d\n', ...
%     STIM.clock, STIM.NIFreq, length(STIM.wave));
ch1=[];

AO.stop 
AO.Rate = STIM.NIFreq; 
AO.TriggersPerRun = 1;
queueOutputData(AO, STIM.wave); % wave is FULL

if isempty(RP)
    RP=actxcontrol('rpco.x', [5 5 26 26]);
    if(RP.ConnectRP2('USB',1) == 0)
        error('Failed to connect to RP2.1');
    end
end
% TagNum = double(RP.GetNumOf('ParTag'));
% for loop = 1:TagNum
%     TagName{loop} = RP.GetNameOf('ParTag', loop)
% end

nRecordPoints = floor(duration*STIM.sample_freq)+1000;
err = rp_setup(nRecordPoints, 'microphone');
startBackground(AO); % get ni board read to go, then trigger the rp
% fprintf(1, "duration: %f  sf: %f  n %d\n", duration, STIM.sample_freq, nRecordPoints);
% return
tic
recdur = nRecordPoints/STIM.sample_freq;
RP.SetTagVal('Rec_Size', nRecordPoints);
RP.SoftTrg(1); 

% Wait until buffer fills

while(toc < recdur+0.1) %
    curindex = RP.GetTagVal('Data_index');
    if curindex >= nRecordPoints-1000
        break
    end
    stf = check_stop(0);
    if stf == 1 % successful STOP from the button
        err = 1;
        return;
    end

end
rp_setup(nRecordPoints, 'stop');
% read the data and stop the hardware
% ch1 = double(RP.ReadTagV('data_out1', 0, nRecordPoints));
% think about this:
%ch2=RP.ReadTagVEX('data_out2', 0, nRecordPoints, 'I32', 'F64', 2);
ch2 = double(RP.ReadTagV('data_out2', 0, nRecordPoints));

AO.stop;
set_attn(-1);
RP.SoftTrg(2);
RP.Halt;
end       