function [Speaker, Mic] = getSpeakerMic()
%getSpeakerMic Read the speaker and mic configuration from the gui
% parse it and return value in global variable.
hs = findobj('tag', 'ABR_Speaker');
ispk = get(hs, 'Value');
switch ispk
    case 1
        Speaker = 'EC1#4101 closed field';
    case 2
        Speaker = 'ES1#3539 free field';
    case 3
        Speaker = 'MF1#1013 free field';
    case 4
        Speaker = 'MF1#1956 closed field';
    case 5
        Speaker = 'FF1#1013 free field (dome)';
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

