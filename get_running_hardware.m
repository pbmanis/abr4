function [HW] = get_running_hardware(HW)
% get the current hardware configuration - are we on a NI/TDT system or not?
% pulled from Acq3.m, 4/30/2010. Paul B. Manis
% Modified 9/26-29 2016 to use new Session mode for NIDAQmx
%
HW.ACQ_DEVICE = 'none'; 
HW.DEVICE_ID = -1; 
notice = '';
switch(computer)
    case {'PC', 'PCWIN', 'PCWIN64'}
        % set up the NI acquisition and hardware interface
        devices = daq.getDevices();
        HW.AO = daq.createSession('ni');
        HW.AdaptorName = 'nidaq';
        HW.ACQ_DEVICE = 'nidaq';
        HW.DEVICE_ID = 1;
        notice = 'NIDAQ found and registered';
%                daqreset % make sure devices are cleared and accessible.
        devchk = addAnalogOutputChannel(HW.AO,'Dev1', 0, 'Voltage');
        addTriggerConnection(HW.AO, 'External', 'Dev1/PFI0', 'StartTrigger');
        HW.AO.Connections(1).TriggerCondition = 'RisingEdge';

        if(~isvalid(devchk))
            HW.ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
            HW.DEVICE_ID = -1;
            notice = 'NIDAQ/WINSOUND not found, running in test mode';
        end
    case {'MAC', 'MACI', 'MACI64'}
        HW.ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
        HW.DEVICE_ID = -1;
        notice = 'Computer is Mac: running in test mode';
    otherwise
        HW.ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
        HW.DEVICE_ID = -1;
        notice = 'Computer type not recognized, running in test mode'; 
end
fprintf(2, '%s\n', notice);

end
     