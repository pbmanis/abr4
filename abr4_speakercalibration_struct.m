classdef abr4_calibration_struct
    % Calibration properties for abr4
    % holds all of the calibration information for abr4
    
    properties
  % calibration structure
        RefSPL
        Freqs
        maxdB
        dBSPL
        dBSPL_bp
        dBSPL_nf
        Vmeas
        Vmeas_bp
        Gain
        CalAttn
        Speaker
        Microphone
        Date
        DateTime
        
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %   
            obj.RefSPL = '';
            obj.Freqs = [];
            obj.maxdB = [];
            obj.dBSPL = [];
            obj.dBSPL_bp = [];
            obj.dBSPL_nf = [];
            obj.Vmeas = 0.;
            obj.Vmeas_bp = 0.;
            obj.Gain = 0.;
            obj.CalAttn = 0.;
            obj.Speaker = '';
            obj.Microphone = '';
            obj.Date = '';
            obj.DateTime = [];
        end
    end
end

