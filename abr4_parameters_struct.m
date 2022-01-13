classdef abr4_parameters_struct
    % plot properties for abr4
    % holds all of the data information for abr4
    
    properties
        % data structure
        ABR4_Pars
        ABR4_FIG
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %
            obj.ABR4_Pars = [];
            obj.ABR4_FIG = [];
            
        end
    end
end

