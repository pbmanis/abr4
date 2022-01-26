function [ch1, ch2, HW, err] = calstim(nRecordPoints, HW, STIM)
% Load and record signal to do a calibration 
% assumes stimulus waveform is already in STIM.wave and that 
% AO and RP are already setup.
% 5/1/2010 Paul B. Manis
% converted to matlab session interface to NIDAQ 9/28/2016
%

HW.AO.stop;
HW.AO.Rate = STIM.NIFreq; 
% queueOutputData(AO, STIM.wave); % wave is FULL
HW.AO.TriggersPerRun = 1;
queueOutputData(HW.AO, STIM.wave); % wave is FULL

if isempty(HW.RP)
    HW.RP=actxcontrol('rpco.x', [5 5 26 26]);
    if(HW.RP.ConnectRP2('USB',1) == 0)
        error('Failed to connect to RP2.1');
    end
end

[HW, err] = rp_setup(HW, STIM, nRecordPoints, 'microphone');
startBackground(HW.AO); % get ni board read to go, then trigger the rp

HW.RP.SoftTrg(1); % start.

TraceDuration = nRecordPoints/STIM.sample_freq;
curindex=HW.RP.GetTagVal('Index');
lastindex = curindex;

pause(TraceDuration)
HW.RP.SoftTrg(2);

% ch1=double(RP.ReadTagV('data_out1', 0, nRecordPoints));
ch1 = 0;
ch2=double(HW.RP.ReadTagV('data_out2', 0, nRecordPoints));
HW.AO.stop;
set_attn(HW, -1);
HW.RP.SoftTrg(2);
HW.RP.Halt;

end       