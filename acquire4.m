function [data, STIM, chdata, err] = acquire4(cmd, HW, STIM, PLOTS, GUI, CALIBRATION, varargin)
% Data acquisition routine (for all acquisition)
% Expects TDT RP2.1 for acquistion, PA5 attenuator, and NI6113 DAC card.
% Includes the ability to perform limited testing using a sound card.

% 7/15/03 p. manis.
% Acquire data through TDT System III with an RP2.1.
% Arms stimulus for fast NI card, and collects navg waveforms
% Requires call to hardware initialization and getstimparams before calling.
%
% 3/29/07 new version.
% samples individual traces, not using avgbuf.
% also samples 2 channels, to monitor microphone output as well.
% modified so that STIM.block allows rapid multiple reps to be collected in
% sequence.
% the number of reps in a block is STIM.np, the interval between reps is
% ipi. the waveform is sent out the 14-bit dac.
%
% 5/25/2007 pbm.
% Modifications to use a different structure for the stimulus parameters.
% now breaks out:
% StimPerSweep, SweepIPI (e.g., for the stimuli that are presented in one
% stimulus sweep)
% NSweeps, InterSweepInterval (e.g., for the sweeps themselves).
% 4-28-2010 P. Manis
% Major changes to improve calibration stuff automatic - 5/2011. P. Manis
%
% 9-27-2016 P. Manis
% Changed to session-based acquisition, per new Matlab style. Now
% works with Matlab 2016.
%
% Jan 13 2021 P. Manis
% Removing globals (new branch 'no-globals'). All data structures are
% defined as 'persistent' in abr4. Structures are defined in classes.
%
%--------------------------------------------------------------------------
%Tessa's variable (saving raw data):
rawdata=[];
rawdatacount=1;
%--------------------------------------------------------------------------
% Calibration Only:
% Check the input cmd argument for commands that invoke calibration-related
% routines.
%--------------------------------------------------------------------------
err = 0;
data = [];
chdata = [];
if nargin > 0 && (strcmp(cmd, 'microphone') ...
        || strcmp(cmd, 'microphone104') ...
        || strcmp(cmd, 'calibrate') ...
        || strcmp(cmd, 'checkcal'))
    if isempty(HW.AO) % only do this if we are using NI
        return;
    end
    correctCal = 0; % do calibration and store to disk.
    if strcmp(cmd, 'checkcal')
        correctCal = 1; % allows us to check the calibration.
        cmd = 'calibrate';
    end
    HW = calibrations(cmd, correctCal, HW, CALIBRATION, STIM, GUI);
    return;
end

if ~isempty(varargin) && strcmp(cmd, 'attn')
    local_attn = varargin{1};
else
    local_attn = 120.0;
end

%--------------------------------------------------------------------------
% Otherwise we do "Regular acquisition":
%--------------------------------------------------------------------------

err = 0;
debugTiming = 0;
NREJECT = 0;

if STIM.StimPerSweep == 1
    nRecordPoints = floor(STIM.avg_dur*STIM.sample_freq/1000.);
else
    nRecordPoints = floor((STIM.StimPerSweep*STIM.ipi/1000)*STIM.sample_freq); % convert from msec to sec, samples/second
end
nStimPoints = floor((STIM.StimPerSweep*STIM.ipi/1000)*STIM.NIFreq);
InterpFreq = 100000;

timebase_Stim = 0:1/STIM.NIFreq:(nStimPoints-1)/STIM.NIFreq; % express rate in msec
timebase_Record = 0:(1/STIM.sample_freq):(nRecordPoints-1)/STIM.sample_freq;
stimBlockLen = floor(STIM.ipi/1000*STIM.NIFreq);
recBlockLen = floor(STIM.ipi/1000.*STIM.sample_freq);
interBlockLen = floor(STIM.ipi/1000.*InterpFreq);
timebase_Interpolate = 0:(1/InterpFreq):((STIM.StimPerSweep*interBlockLen-1)/InterpFreq);
timebase_Display = 1000*(0:(1/InterpFreq):(interBlockLen-1)/InterpFreq);

replot = 0;
t0 = 1;
t1 = floor(1e-3*0.5*STIM.sample_freq);
t2 = 1;
t3 = floor(1e-3*2.0*STIM.sample_freq);
tbs = 1;
tbe = floor(1e-3*0.9*STIM.sample_freq);

tp0 = floor(1e-3*2.0*STIM.sample_freq);
tp1 = floor(1e-3*6.0*STIM.sample_freq);



% although we generated a full stimulus (or used to), we're just going to get the first
% stimulus waveform, or the first 2 if alternating, and let the dac run out
% repeatedly.
if(STIM.Alternate == 0)
    stdur = STIM.ipi;
    stdacrepeat = STIM.StimPerSweep;
    stpmax = floor((stdur*10^-3)/(1/STIM.NIFreq));
    aosamples = length(STIM.wave);
else
    stdur = 2.0*STIM.ipi;
    stdacrepeat = ceil(STIM.StimPerSweep/2);
    stpmax = floor((stdur*10^-3)/(1/STIM.NIFreq));
    aosamples = length(STIM.wave);
end
if ~isempty(HW.AO) % only do this if we are using NI
    if(aosamples < stpmax)
        stpmax = aosamples;
    end
    swave = STIM.wave(1:stpmax);
else
    swave = STIM.wave;
end
scaleFactor_ch1 = 1.0/STIM.amp_gain; % (10/(2^23))*1000000/STIM.amp_gain; % convert to microvolts, correcting for a/d in rp2.1
scaleFactor_ch2 = 1.0; % (10/(2^23)); % convert to millivolts, correcting for a/d in rp2.1

skip = 1;

hrej = findobj('tag', 'ABR_Nreject');
if(~isempty(hrej))
    set(hrej, 'string', sprintf('%d', NREJECT));
end
record_dur = max(timebase_Record);
stim_dur = max(timebase_Stim);

if isempty(HW.RP)
    ar = audiorecorder(44100, 16, 1); % use audio channel for input
end
if isempty(HW.AO)
    PL = audioplayer(swave, 44100);
end

data = zeros(4, interBlockLen);  % cumulative data across sweeps, neg, pos, ch1 and ch2.
Total_np = 0;
Total_nn = 0;
Total_all = 0;

% this routine accumulates the sum of the good trials into data
% we average the data at the end, using nn and np (the number of good + and -
% trials)
%data(1,:) main signal on in-1,
%data(2,:) alternate polarity
%data(3,:)  in-2 - audio signal
%data(4,:)  in-2 - audio signal
plotSweep = cell(3, STIM.NSweeps);
set_status('Running');
for i_sweep = 1:STIM.NSweeps % loop over all the sweeps.
    % check for the stop button
    state = check_status(GUI);
    if strcmp(state, 'Stopped')
        abr4('Stop', 'button in sweep loop');
        err = 1;
        return
    end
    sweepdata = zeros(4, interBlockLen);
    Sweep_np = 0;
    Sweep_nn = 0;
    Sweep_nonalt = 0;
    
    %    fprintf(2,'sweep: %d   stdacrepeat: %d \n', i, stdacrepeat);
    STIM.Monitor = get_SignalPlotFlag(2); % allow plot flag to change when running
    HW.IN_ACQ = true; %#ok<*NASGU>
    sweepreject = 0;
    tic % time mark
    HW = set_attn(HW, local_attn);
    pause(0.01);
    HW.AO.stop
    HW.AO.Rate = STIM.NIFreq;
    queueOutputData(HW.AO, STIM.wave); % wave is FULL
    HW.AO.TriggersPerRun = 1; % set(AO, 'repeatoutput',  1);
    set_status('Running');
    tic
    startBackground(HW.AO); % get ni board read to go, then trigger the rp
    pause(0.1);  % give the system time to arm
    [HW, err] = rp_setup(HW, STIM, nRecordPoints+1000, 'Start');
    if err == 1
        fprintf(2, 'Hardware failed to start\n');
        set_status('Stopped');
        abr4('Stop', 'Error in RP setup');
        return;
    end
    recdur = nRecordPoints/STIM.sample_freq;
    curindex=HW.RP.GetTagVal('Data_index');
    % Wait until buffer fills
    while( toc < recdur) %
        curindex=HW.RP.GetTagVal('Data_index');
        pause(0.01);
        state = check_status(GUI);
        if strcmp(state, 'Stopped')
            abr4('Stop', 'button in toc loop');
            err = 1;
            return
        end
        % read the data and stop the hardware
        pause(1.5); % wait for the data to become fully available
        ch1 = double(HW.RP.ReadTagV('data_out1', 0, nRecordPoints));
        % think about this:
        %        da=HW.RP.ReadTagVEX('data_out1', 0, nRecordPoints, 'I32', 'F64', 2);
        if STIM.Monitor
            ch2 = double(HW.RP.ReadTagV('data_out2', 0, nRecordPoints));
        end
        state = check_status(GUI);
        if strcmp(state, 'Stopped')
            abr4('Stop', 'button in toc loop');
            err = 1;
            return
        end
    end
    [HW, ~] = rp_setup(HW, STIM, nRecordPoints, 'Stop');

    set_status('Waiting');
    ch1=ch1*scaleFactor_ch1; % express in microvolts.
    %%%%%%%%%%%%%%%%%%%%%%
    %collecting the raw data
%     rawdata(rawdatacount,:)=ch1;
%     rawdatacount=rawdatacount+1;
    %%%%%%%%%%%%%%%%%%%%%%
    
    % now split according to the number of stimulus pulses in the waveform
    ch1i = interp1(timebase_Record, ch1, timebase_Interpolate, 'linear', 'extrap');
    ch1r = reshape(ch1i(1:(interBlockLen*STIM.StimPerSweep)), ...
        [interBlockLen, STIM.StimPerSweep]);
    
    c=clock;
    fname = sprintf('%4d%02d%02d-%02d%02d-raw-%d.txt', c(1), c(2), c(3), c(4), c(5), i_sweep);
    % save(fname, 'ch1r', '-ascii', '-tabs');
    %   fprintf(2, 'max cha: %f nV\n', max(abs(ch1*1e9)));
    
    if(STIM.Monitor == 1)
        ch2 = ch2*scaleFactor_ch2; % express in volts.
        ch2i = interp1(timebase_Record, ch2, timebase_Interpolate, 'linear', 'extrap');
        ch2r = reshape(ch2i(1:(interBlockLen*STIM.StimPerSweep)), ...
                [interBlockLen, STIM.StimPerSweep]);
            ch2m = mean(ch2r, 2);
    else
        ch2 = [];
    end
    %----------------------------------------------------------------------
    % Artifact rejection
    % Reject traces on two potential causes:
    % 1. The baseline is noisy.
    % 2. the absolute value of the signal is greater than
    %    STIM.avg_reject, in microvolts.
    % The array iy holds the indices into the samples in the sweep
    % that meet these criteria.
    %----------------------------------------------------------------------
    rejected_indices = [];
    reject_on_amplitude = [];
    reject_on_std = [];
    response_indices = 1:STIM.StimPerSweep;
    if STIM.avg_reject > 0
        [~, reject_on_amplitude] = find(abs(ch1r) > STIM.avg_reject*10^-6);
        allstds=std(ch1r(tbs:tbe,:));
        w = mean(allstds);
        reject_on_std = find(allstds > w*STIM.rms_reject)';
        rejected_indices = [reject_on_amplitude; reject_on_std]; % concatenate...
        rejected_indices = unique(rejected_indices);  % only count each once
    end

    response_indices = setxor(response_indices, rejected_indices); % exclude those that were rejected.
    sweeps_rejected = length(rejected_indices);
    NREJECT = NREJECT + sweeps_rejected;
    hrej = findobj('tag', 'ABR_Nreject');
    if(~isempty(hrej))
        set(hrej, 'string', sprintf('N=%d [SW: %d]', NREJECT, sweeps_rejected));
    end

    for ii = 1:length(response_indices) % break out the data from the acquisition to individual arrays
        if find(rejected_indices == response_indices(ii))  % If response was rejected, leave it out
            continue
        end
        if(STIM.Alternate)
            switch(mod(ii,2)) % put alternate data into different data sets.
                case 1
                    sweepdata(1,:) = sweepdata(1,:) + ch1r(:, response_indices(ii))'; % sum alternate polarites into different channels
                    data(1,:)  = data(1,:) + ch1r(:, response_indices(ii))';
                    Sweep_np = Sweep_np + 1;
                    Total_np = Total_np + 1;
                    if(STIM.Monitor)
                        sweepdata(3,:) = sweepdata(3,:) + ch2r(:, response_indices(ii))';
                        data(3,:) = data(3,:) + ch2r(:, response_indices(ii))';
                    end % input channel 2 is microphone monitor of STIM...
                case 0
                    sweepdata(2,:) = sweepdata(2,:) + ch1r(:, response_indices(ii))'; % the alternation group
                    data(2,:)  = data(2,:) + ch1r(:, response_indices(ii))';
                    Sweep_nn = Sweep_nn + 1;
                    Total_nn = Total_nn + 1;
                    if(STIM.Monitor)
                        sweepdata(4,:) = sweepdata(4,:) + ch2r(:, response_indices(ii))';
                    end % input channel 2 is microphone monitor of STIM...
                otherwise
            end
        else % not in alternation mode
            sweepdata(1,:) = sweepdata(1,:) +  ch1r(:, response_indices(ii))';
            data(1,:) = data(1,:) + ch1r(:, response_indices(ii))';
            if(STIM.Monitor)
                sweepdata(3,:) = sweepdata(3,:) + ch2r(:, response_indices(ii))';
            end
            Sweep_nonalt = Sweep_nonalt + 1;
            Total_all = Total_all + 1;
        end
    end

    set(GUI.hrep, 'string', sprintf('%d', i_sweep*STIM.StimPerSweep)); % update acquisition trial counter
    if(mod(replot, 10) == 0)
        replot = 0;
    else
        replot = replot + 1;
    end
    % average the data from accepted sweeps, then sum in to the data arrays
    if(STIM.Alternate)
        sweepdata(1,:) = sweepdata(1,:)/Sweep_np;
        sweepdata(2,:) = sweepdata(2,:)/Sweep_nn;
        if STIM.Monitor
            sweepdata(3,:) = sweepdata(3,:)/Sweep_np;
            sweepdata(4,:) = sweepdata(4,:)/Sweep_nn;
            data(3,:) = data(3,:) + sweepdata(3,:);
            data(4,:) = data(4,:) + sweepdata(4,:);
        end
        set(PLOTS.data, 'clipping','On');
        hold on; % set(PLOT.PLOTS.data, 'hold', 'On');
        runningAvg = ((data(1,:)/Total_np)+(data(2,:)/Total_nn))/2.0;
        runningAvg = runningAvg - mean(runningAvg(tbs:tbe));
        sweepAvg =(sweepdata(1,:)+sweepdata(2,:))/2.0;
        sweepAvg = sweepAvg - mean(sweepAvg(tbs:tbe));

        plotSweep{1,i_sweep} = plot(PLOTS.data, timebase_Display, sweepAvg);
        plotSweep{1,i_sweep}.Color=[0.6, 0.3, 0.3, 0.5];
        
        if i_sweep == 1
            avgplot = plot(PLOTS.data, timebase_Display, runningAvg, 'k-', 'LineWidth', 2.0);
        else
            set(avgplot,'Xdata',timebase_Display,'YData',runningAvg);
        end
        md = max(abs(sweepdata(1,:)));
        
    else
        sweepdata(1,:) = sweepdata(1,:)/Sweep_nonalt;
        sweepdata(1,:) = sweepdata(1,:) - mean(sweepdata(1,tbs:tbe));
        runningAvg = data(1,:)/Total_all;
        runningAvg = runningAvg - mean(runningAvg(1, tbs:tbe));
        if STIM.Monitor
            sweepdata(3,:) = sweepdata(3,:)/Sweep_nonalt;
            data(3,:) = data(3,:) + sweepdata(3,:);
        end
        
        plotSweep{1, i_sweep} = plot(PLOTS.data, timebase_Display, sweepdata(1,:), 'b-');
        if i_sweep == 1
            avgplot = plot(PLOTS.data, timebase_Display, runningAvg, 'k-', 'LineWidth', 2.0);
        else
            set(avgplot,'Xdata',timebase_Display,'YData',runningAvg);
        end
        plotSweep{1,i_sweep}.Color=[0.6, 0.6, 0.6, 0.5];
        md = max(abs(sweepdata(1,:)));
    end
  
    set_scale(PLOTS)

    set(PLOTS.data, 'xlim', [0 STIM.avg_dur]);
    set(PLOTS.data, 'clipping', 'on');
    drawnow;
    
    if(STIM.Monitor == 1) % only do this if we have the monitor set on
        hold off
        plot(PLOTS.signal2, timebase_Display, sweepdata(3,:), 'color', 'red');
        if(STIM.Alternate)
            hold on
            plot(PLOTS.signal2, timebase_Display, sweepdata(4,:), 'color', 'cyan');
        end
        md = max(abs(sweepdata(3,:)));
        set(PLOTS.signal2, 'xLimMode', 'Auto');
        set(PLOTS.signal2, 'yLimMode', 'Auto');
    end
    
    % delay to next sweep:
    if i_sweep < STIM.NSweeps % (but not on last sweep)
%        fprintf(2, 'i: %d n: %d\n', i_sweep, STIM.NSweeps);
        while toc < STIM.InterSweepInterval
            pause(0.01);
            state = check_status(GUI);
            if state == "Stopped"
                abr4('Stop', 'in delay to next sweep');
                err = 1;
                return
            end
        end
    end
end  % END of the Sweeps FOR loop

% divde accumulated averages by the number of sweeps accepted.
% subtract the baseline

if STIM.Alternate
    data(1,:) = data(1,:)/Total_np;
    data(1,:) = data(1,:) - mean(data(1, tbs:tbe));
    data(2,:) = data(2,:)/Total_nn;
    data(2,:) = data(2,:) - mean(data(2, tbs:tbe));
else
    data(1,:) = data(1,:)/Total_all;
    data(1,:) = data(1,:) - mean(data(1, tbs:tbe));
    data(2,:) = data(2,:)/Total_all; 
    data(2,:) = data(2,:) - mean(data(2, tbs:tbe));
end

if STIM.Monitor
    data(3,:) = data(3,:)/STIM.NSweeps;
    data(4,:) = data(4,:)/STIM.NSweeps;
end
chdata = [ch1, ch2];

STIM.ACQPars_np = Total_np;
STIM.ACQPars_nn = Total_nn;
STIM.ACQPars_nall = Total_all;
STIM.ACQPars_tb = timebase_Display;
set_status('Stopped');  % was 'Done'
% fid = fopen('Tessaraw.txt','w');
% fprintf(fid,'%6.2f  %12.8f\n',rawdata);
% fprintf(2, 'max ch 2: %f\n', md);

end


