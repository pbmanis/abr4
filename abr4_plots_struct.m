classdef abr4_plots_struct
    % plot properties for abr4
    % holds all of the data information for abr4
    
    properties
        % data structure
        PLOTHANDLES
        responsemap
        quiver
        data
        signal1
        signal2
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instance of this class
            %   
            obj.PLOTHANDLES = [];
            obj.responsemap = [];
            obj.quiver = [];
            obj.data = [];
            obj.signal1 = [];
            obj.signal2 = [];
        end
    end
end

