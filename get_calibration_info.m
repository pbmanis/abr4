function [CALIBRATION] = get_calibration_info(Speaker, Mic, CALIBRATION)

% verify dates of last calibrations. If more than one week, refuse to
% proceed.
mic_struct = abr4_microphone_struct();
mic = load(sprintf('microphone_%s.cal', Mic), '-mat');
MIC = mic_struct.from_struct(mic.MIC);

spkr = load(sprintf('frequency_%s.cal', Speaker), '-mat');
SPKR = spkr.CAL;
chk = load(sprintf('chk75db_%s.cal', Speaker), '-mat');
CHK75 = chk.CHKCAL;
fprintf(1, "Using Calibrations from:\n");
fprintf(1, "   Microphone:   %24s  on %s\n", Mic, MIC.Date);
fprintf(1, "   Speaker Cal:  %24s  on %s\n", Speaker, SPKR.Date);
fprintf(1, "   Speaker 75dB: %24s  on %s\n", Speaker, CHK75.Date);
today = datetime;
mic_date = datetime(MIC.Date);
spk_date = datetime(SPKR.Date);
chk_date = datetime(CHK75.Date);
mic_days = caldays(between(mic_date, today, 'days'));
CALIBRATION.Needs_Cal = false;
if mic_days > 7
    fprintf(2, 'Last microphone calibration was more than 1 week ago, please recalibrate first\n');
    CALIBRATION.Needs_Cal = true;
end
spk_days = caldays(between(spk_date, today, 'days'));
if spk_days > 7
    fprintf(2, 'Last speaker calibration was more than 1 week ago, please recalibrate first\n');
    CALIBRATION.Needs_Cal = true;
end
chk_days = caldays(between(chk_date, today, 'days'));
if chk_days > 7
    fprintf(2, 'Last check75dB calibration was more than 1 week ago, please recalibrate first\n');
    CALIBRATION.Needs_Cal = true;
end
if mic_days ~= spk_days || spk_days ~= chk_days
    fprintf(2, 'Calibrations were not done on same day; please recalibrate first\n');
else
    fprintf(1, 'All calibrations are from the same day, OK to proceed\n');
end
CALIBRATION.MIC = MIC;
CALIBRATION.SPKR = SPKR;
CALIBRATION.CHK75 = CHK75;

end