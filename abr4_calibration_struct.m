classdef abr4_calibration_struct
    % Calibration properties for abr4
    % holds all of the calibration information for abr4
    
    properties
  % calibration structure
       MIC  % holds microphone struct from calibrations being used in run
       CHK75  % holds chk75db structure from calibratoin
       SPKR  % holds speaker calibration. 
       Needs_Cal
       SPLCAL
       
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %   
            obj.MIC = [];
            obj.CHK75 = [];
            obj.SPKR = [];
            obj.Needs_Cal = false;
            obj.SPLCAL = [];
        end
    end
end

