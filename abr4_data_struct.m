classdef abr4_data_struct
    % data properties for abr4
    % holds all of the data information for abr4
    
    properties
        % data structure
        DATAa = [];
        DATAp = [];
        DATAn = [];
        CDATAp = [];
        CDATAn = [];
        CHDATA = [];
        C_CHDATA = []
        Timebase = [];
        REFERENCE = [];
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %   
            obj.DATAa = [];  % all
            obj.DATAp = [];
            obj.DATAn = [];
            obj.CDATAp = [];
            obj.CDATAn = [];
            obj.CHDATA = [];
            obj.C_CHDATA = [];  % Raw data
            obj.REFERENCE = [];
            obj.Timebase = [];
        end
    end
end

