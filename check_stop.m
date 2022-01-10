function [stopindicator] = check_stop(stopflag)

global AO RP STOP IN_ACQ

stopindicator = 0;

if stopflag == 1 % force stop
    STOP = 1;
end;

if(STOP == 1) % while waiting, check for stop.
    RP.SoftTrg(0); % invoke(RP, 'softtrg', 0);
    RP.Halt; % invoke(RP, 'halt');
    set_attn(120);
    AO.stop;
    %queueOutputData(AO, 0); % make sure output is really back to 0
    IN_ACQ = 0;
    stopindicator = 1; % signal that we stopped
    STOP = 0;
    set_status('Stopped');
    return;
end
end
