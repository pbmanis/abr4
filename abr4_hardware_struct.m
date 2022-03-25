classdef abr4_hardware_struct
    % hardware properties for abr4
    % holds all of the hardware structure information for abr4
    
    properties
  % Hardware structure
        AdaptorName
        HARDWARE
        DEVICE_ID
        ACQ_DEVICE
        AO
        RP
        PA5
        STOP_BUTTON_HIT
        IN_ACQ
        STIM
    end
    
    methods
        function obj = initialize(obj)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.AdaptorName = '';
            obj.HARDWARE = '';
            obj.DEVICE_ID = [];
            obj.ACQ_DEVICE = [];
            obj.AO = [];
            obj.RP = [];
            obj.PA5 = [];
            obj.STOP_BUTTON_HIT = false;  % Set true when stop button is hit, set false when stopped
            obj.IN_ACQ = false;  % True when acquisition is up an running
            obj.STIM = [];
            
        end
    end
end

