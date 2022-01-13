function [Speaker, Mic] = getSpeakerMic()
%getSpeakerMic Read the speaker and mic configuration from the gui
% parse it and return value in global variable.
hs = findobj('tag', 'ABR_Speaker');
ispk = get(hs, 'Value');
ispk = ispk(1);

switch ispk
    case 1
        Speaker = 'EC1';
    case 2
        Speaker = 'ES1';
    case 3
        Speaker = 'MF1';
    otherwise
        Speaker = 'Unknown';
end
hs = findobj('tag', 'ABR_Microphone');
imic = get(hs, 'Value');

switch imic
    case 1
        Mic = '7012#39279'; % 1/2 microphone
    case 2
        Mic = '7016#9945'; % 1/4 microphone
    case 3
        Mic = '7016#10252'; % 1/4 microphone (2018)
    otherwise
        Mic = 'Unknown';
end

end

