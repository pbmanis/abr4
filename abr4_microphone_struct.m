classdef abr4_microphone_struct
    % data properties for abr4
    % holds all of the data information for abr4
    
    properties
        % data structure
        Vrms 
        Vref_c  % convert to RMS
        Vref_bp
        Microphone
        Date
        dBPerVPa
        mVPerPa
        
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %   
            obj.Vrms = [];
            obj.Vref_c = [];
            obj.Vref_bp = [];
            obj.Microphone = '';
            obj.Date = '';
            obj.dBPerVPa = [];
            obj.mVPerPa = [];
        end
        
        function obj = from_struct(obj, s)
            obj.Vrms = s.Vrms;
            obj.Vref_c = s.Vref_c;
            obj.Vref_bp = s.Vref_bp;
            obj.Microphone = s.Microphone;
            obj.Date = s.Date;
            obj.dBPerVPa = s.dBPerVPa;
            obj.mVPerPa = s.mVPerPa;
        end
    end
end

