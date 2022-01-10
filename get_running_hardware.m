function get_running_hardware()
%% get the current hardware configuration - are we on a NI/TDT system or not?
%% pulled from Acq3.m, 4/30/2010. Paul B. Manis
%% Modified 9/26-29 2016 to use new Session mode for NIDAQmx
%%

global AO DEVICE_ID ACQ_DEVICE
ACQ_DEVICE = 'none'; %#ok<NASGU>
DEVICE_ID = -1; %#ok<NASGU>
notice = ''; %#ok<NASGU>
switch(computer)
    case {'PC', 'PCWIN', 'PCWIN64'}
        % set up the NI acquisition and hardware interface
        devices = daq.getDevices();
        AO = daq.createSession('ni');
        hw.AdaptorName = 'nidaq';
        ACQ_DEVICE = 'nidaq';
        DEVICE_ID = 1;
        notice = 'NIDAQ found and registered';
%                daqreset % make sure devices are cleared and accessible.
        devchk = addAnalogOutputChannel(AO,'Dev1', 0, 'Voltage');
        addTriggerConnection(AO, 'External', 'Dev1/PFI0', 'StartTrigger');
        AO.Connections(1).TriggerCondition = 'RisingEdge';

        if(~isvalid(devchk))
            ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
            DEVICE_ID = -1;
            notice = 'NIDAQ/WINSOUND not found, running in test mode';
        end
    case {'MAC', 'MACI', 'MACI64'}
        ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
        DEVICE_ID = -1;
        notice = 'Computer is Mac: running in test mode';
    otherwise
        ACQ_DEVICE = 'none'; % no hardware - we have a dummy acquisition mode
        DEVICE_ID = -1;
        notice = 'Computer type not recognized, running in test mode'; %#ok<NASGU>
end;
fprintf(2, '%s\n', notice);

end
     