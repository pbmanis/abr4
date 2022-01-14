classdef abr4_hardware_struct
    % hardware properties for abr4
    % holds all of the hardware structure information for abr4
    
    properties
  % Hardware structure
        HARDWARE
        DEVICE_ID
        ACQ_DEVICE
        AO
        RP
        PA5
        STOP
        IN_ACQ
        STIM
    end
    
    methods
        function obj = initialize(obj)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.HARDWARE = '';
            obj.DEVICE_ID = [];
            obj.ACQ_DEVICE = [];
            obj.AO = [];
            obj.RP = [];
            obj.PA5 = [];
            obj.STOP = false;
            obj.IN_ACQ = false;
            obj.STIM = [];
            
        end
    end
end

