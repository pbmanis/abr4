function [ err ] = verifyConfiguration(mode, cmd)
% Manual verification of configuration for calibration
% or for recording setup.

err = 0;
msg = {};
msg{1} = sprintf('Please check the following: \n');
n = 2;
[Speaker, Mic] = getSpeakerMic();
fprintf(2, "Mic: %s \n", Mic);
switch mode
    case 'calibrate'
        msg{n} = sprintf('Microphone: \n  Microphone type = %s', Mic);
        n = n + 1;
        msg{n} = sprintf('   Microphone Preamp Gain is 20dB');
        n = n + 1;
        switch cmd
            case "microphone"
                msg{n} = sprintf('   Calibrator Standard is set to 94dBSPL and ON');
                n = n + 1;
            case 'microphone104'
                msg{n} = sprintf('   Calibrator Standard is set to 104dBSPL and ON');
                n = n + 1;
            otherwise
        end
    otherwise
end
msg{n} = sprintf('\nSpeaker: \n   Speaker is %s', Speaker);
n = n + 1;
msg{n} = sprintf('   Speaker driver attenuation is -6dB');
n = n + 1;
msg{n} = sprintf('\nAre your settings Correct?');
button = questdlg(msg,'Hardware Settings','Yes','No','Yes') ;
if strcmp(button,'No')
    err = 1;
end

end

