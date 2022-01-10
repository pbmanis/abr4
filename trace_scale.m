function trace_scale(md, hdata)
%trace_scale provides logical "autoscaling" with steps
%   
    maxd = [1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, ...
            1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, ...
            1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, ...
            1e-2, 5e-2, 1e-1, 5e-1, 1, 2, 5, 10];
    mind = - maxd;
    imd = find(maxd > md);
    imind = find(mind < md);
    
    if (~isempty(imd))
        set(hdata, 'YLim', [-maxd(imd(1)) maxd(imd(1))]);
    else
    %    set(hdata, 'YLimMode', 'auto');
    end;

end

