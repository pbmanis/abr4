classdef abr4_gui_struct
    % GUI handles for abr4
    % holds all of the GUI handles for abr4
    
    properties
       hstat
       hsr  % handle to graphics
        hchk
        hisi  % handle to graphics
        hrep  % handle to graphics
        hsps
        hipi
        hdur
        hreject
        hreject2
        hgain
        hmin
        hmax
        hstep
        hflist
        hcursor
        hstimfilename
        hcurrentfrequency
        hcurrentSPL
        hToneMapDuration
    end
    
    methods
        function obj = initialize(obj)
            % Construct an instfance of this class
            %
            obj.hstat = findobj('Tag', 'ABR_Status');
            obj.hsr = findobj('Tag', 'ABR_NSweeps'); % number of sweeps to average across
            obj.hchk = findobj('Tag', 'ABR_AlternatePolarity');
            obj.hrep = findobj('Tag', 'ABR_StimRep');% rep counter
            obj.hisi = findobj('Tag', 'ABR_InterSweepInterval'); % number of sweeps to average across
            obj.hsps = findobj('Tag', 'ABR_StimPerSweep'); % read interstimulus interval within a sweep
            obj.hipi = findobj('Tag', 'ABR_SweepIPI'); % read interstimulus interval within a sweep
            obj.hdur = findobj('Tag', 'ABR_Duration');% Display duration
            obj.hreject = findobj('Tag', 'ABR_Reject');% artificat rejection criterion: AMplitude criteria
            obj.hreject2 = findobj('Tag', 'ABR_Reject2'); % artificat rejection criterion (stdev criteria)
            obj.hgain = findobj('Tag', 'ABR_Gain');% amplifier gain setting - nominal is 10000.
            obj.hmin = findobj('Tag', 'ABR_MinAttn');
            obj.hmax = findobj('Tag', 'ABR_MaxAttn');
            obj.hstep = findobj('Tag', 'ABR_AttnStep');
            obj.hflist = findobj('Tag', 'ABR_FreqList');
            obj.hcursor = findobj('Tag', 'ABR_cursor');
            obj.hstimfilename = findobj('Tag', 'ABR_StimFile');
            obj.hcurrentfrequency = findobj('Tag', 'ABR_CurrFreq');
            obj.hcurrentSPL = findobj('Tag', 'ABR_CurrSPL');
            obj.hToneMapDuration = findobj('Tag', 'ABR_ToneMapDur');

        end
        function print_handles(obj)
        fprintf(1, "GUI Handles\n");
        disp(obj.hflist)
        fllist = get(obj.hflist, 'string');
        fprintf(1, "fllist: %s", fllist);
        disp(obj.hsr)
        hsr = get(obj.hsr, 'string');
        disp(hsr)
        fprintf(1, "hsr: %s", hsr);
        end
        function [x] = get_handle(obj)
        x = findobj('Tag', 'ABR_FreqList');
        disp('get tag for freqlsit')
        disp(x)
        end
    end
end

