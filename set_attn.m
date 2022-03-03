function [HW] = set_attn(HW, attn)
% set_attn sets the PA5 Programmable attenuator and displays the setting.
% If this is called and the HW.PA5 control has not yet been initialized,
% it will take care of that initialization.
if(attn > 120.0 || attn < 0.0) % out of range : set to maximum attenuation
    attn = 120.0;
end

usbchan = 1;
if(isempty(HW.PA5) || ~ishandle(HW.PA5))
    fprintf(2, 'Attempting to connect to PA5 Attenuators: ');
    HW.PA5=actxcontrol('PA5.x', [5 5 26 26]);
    if(HW.PA5.ConnectPA5('USB', usbchan) == 0)
        fprintf(2, ' ... failed to connect to PA5 %d\n', usbchan);
        return;
    end
    fprintf(2, 'PA5 Connection OK\n');
end
HW.PA5.SetAtten(attn);
hattn = findobj('Tag', 'ABR_CurrAttn');
if(~isempty(hattn))
    set(hattn, 'String', sprintf('%5.1f dB', attn));
end
return;