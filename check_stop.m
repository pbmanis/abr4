function [HW, stopindicator] = check_stop(HW, stopflag)
stopindicator = 0;

if stopflag == 1 % force stop
    HW.STOP = 1;
end

if(HW.STOP == 1) % while waiting, check for stop.
    HW.RP.SoftTrg(0); % invoke(RP, 'softtrg', 0);
    HW.RP.Halt; % invoke(RP, 'halt');
    set_attn(HW, 120);
    HW.AO.stop;
    %queueOutputData(AO, 0); % make sure output is really back to 0
    HW.IN_ACQ = 0;
    stopindicator = 1; % signal that we stopped
    HW.STOP = 0;
    set_status('Stopped');
    return;
end
end
