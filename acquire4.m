function [data, STIM, err] = acquire4(cmd, HW, STIM, PLOTS, GUI, CALIBRATION, varargin)
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
    HW = calibrations(cmd, correctCal, HW, CALIBRATION, STIM);
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

data = zeros(4, interBlockLen);
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

set_status('Running');
plotSweep = cell(3, STIM.NSweeps);
% STIM.NSweeps = 2;
for i = 1:STIM.NSweeps % loop over all the sweeps.
    sweepdata = zeros(4, interBlockLen);
    %    fprintf(2,'sweep: %d   stdacrepeat: %d \n', i, stdacrepeat);
    STIM.Monitor = get_SignalPlotFlag(2); % allow plot flag to change when running
    IN_ACQ = 1; %#ok<*NASGU>
    STOP = 0; % make sure stop is reset first
    sweepreject = 0;
    tic % time mark
    if ~isempty(HW.AO)
        [HW, stf] = check_stop(HW, 0);
        if stf == 1 % successful STOP from the button
            fprintf(2, "stop with button");
            err = 1;
            return;
        end
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
        [HW, err] = rp_setup(HW, STIM, nRecordPoints+1000, 'start');
        if err == 1
            fprintf(2, 'Hardware failed to start\n');
            return;
        end
        recdur = nRecordPoints/STIM.sample_freq;
        curindex=HW.RP.GetTagVal('Data_index');
        % Wait until buffer fills
        while( toc < recdur) %
            curindex=HW.RP.GetTagVal('Data_index');
            pause(0.01);
            [HW, stf] = check_stop(HW, 0);
            if stf == 1 % successful STOP from the button
                fprintf(2, "stop with button");
                err = 1;
                return;
            end
        end
        % read the data and stop the hardware
        pause(0.05); % wait for the data to become fully available
        ch1 = double(HW.RP.ReadTagV('data_out1', 0, nRecordPoints));
        % think about this:
        %        da=HW.RP.ReadTagVEX('data_out1', 0, nRecordPoints, 'I32', 'F64', 2);
        if STIM.Monitor
            ch2 = double(HW.RP.ReadTagV('data_out2', 0, nRecordPoints));
        end
        [HW, ~] = rp_setup(HW, STIM, nRecordPoints, 'stop');
%         fprintf(2, "rpsetup return code: %d", rp_return);
    else
        play(PL); % start it...
        while(~strcmp(PL.running, 'off'))
            [HW, stf] = check_stop(HW, 0);
            if stf == 1
                PL.stop;
                HW.IN_ACQ = false;
                err = 1;
                set_status('Aborted');
                return;
            end
        end
        recordblocking(ar, nRecordPoints/STIM.sample_freq)
        ch1 = getaudiodata(ar);
    end
    HW.IN_ACQ = false;
    ch1=ch1*scaleFactor_ch1; % express in microvolts.
    %%%%%%%%%%%%%%%%%%%%%%
    %collecting the raw data
    rawdata(rawdatacount,:)=ch1;
    rawdatacount=rawdatacount+1;
    %%%%%%%%%%%%%%%%%%%%%%
    
    % now split according to the number of stimulus pulses in the waveform
    ch1i = interp1(timebase_Record, ch1, timebase_Interpolate, 'linear', 'extrap');
    ch1r = reshape(ch1i(1:(interBlockLen*STIM.StimPerSweep)), ...
        [interBlockLen, STIM.StimPerSweep]);
    
    %{
    if(debugTiming) % this is very useful for debugging the timing...
        hf = findobj('Tag', 'debugfig');
        if isempty(hf)
            hf = figure;
            set(hf, 'Tag', 'debugfig');
        end
        figure(hf);
        subplot(2,1,1);
        hold on;
        plot(ch1, 'g');
    end
    %}
        
    c=clock;
    fname = sprintf('%4d%02d%02d-%02d%02d-raw-%d.txt', c(1), c(2), c(3), c(4), c(5), i);
    % save(fname, 'ch1r', '-ascii', '-tabs');
    %   fprintf(2, 'max cha: %f nV\n', max(abs(ch1*1e9)));
    
    if(STIM.Monitor == 1)
        fprintf(1, 'CH2 max in loop: %f, scalefactor: %f\n', max(ch2), scaleFactor_ch2);
        ch2 = ch2*scaleFactor_ch2; % express in volts.
        ch2i = interp1(timebase_Record, ch2, timebase_Interpolate, 'linear', 'extrap');
        ch2r = reshape(ch2i(1:(interBlockLen*STIM.StimPerSweep)), ...
            [interBlockLen, STIM.StimPerSweep]);
        ch2m = mean(ch2r, 2);
        % fprintf(1, 'CH2 scaled and max: %f\n', max(ch2r));
        %         if(debugTiming) % this is very useful for debugging the timing...
        %             hf = findobj('Tag', 'debugfig');
        %             if isempty(hf)
        %                 hf = figure;
        %                 set(hf, 'Tag', 'debugfig');
        %             end;
        %             figure(hf);
        %             subplot(2,1,2);
        %             hold on;
        %             plot(ch2, 'r');
        %         end;
    end
    
    %----------------------------------------------------------------------
    % Artifact rejection code:
    %----------------------------------------------------------------------
    iy = [];
    li = 1:STIM.StimPerSweep;
    
    if STIM.avg_reject > 0
        [~, iy] = find(abs(ch1r) > STIM.avg_reject*10^-6);
        allstds=std(ch1r(tbs:tbe,:));
        w=mean(allstds);
        iy2 = find(allstds > w*STIM.rms_reject)';
        iy = [iy; iy2]; % concatenate...
        iy = unique(iy);
    end
    lin = setxor(li, iy); % exclude those that were rejected.
    NREJECT = NREJECT + length(iy);
    sweepreject = length(iy);
    hrej = findobj('tag', 'ABR_Nreject');
    if(~isempty(hrej))
        set(hrej, 'string', sprintf('%d [%d]', NREJECT, sweepreject));
    end
    Sweep_np = 0;
    Sweep_nn = 0;
    Sweep_nonalt = 0;
    for ii = 1:length(lin) % break out the data from the acquisition to individual arrays
        if(STIM.Alternate)
            switch(mod(ii,2)) % put alternate data into different data sets.
                case 1
                    % fprintf(1, 'Alternate into data 1 %d\n', ii);
                    sweepdata(1,:) = sweepdata(1,:) + ch1r(:, ii)'; % sum alternate polarites into different channels
                    if(STIM.Monitor)
                        sweepdata(3,:) = sweepdata(3,:) + ch2r(:, ii)';
                    end % input channel 2 is microphone monitor of STIM...
                    Sweep_np = Sweep_np + 1;
                case 0
                    sweepdata(2,:) = sweepdata(2,:) + ch1r(:, ii)'; % the alternation group
                    % fprintf(1, 'Alternate into data 2 %d\n', ii);
                    if(STIM.Monitor)
                        sweepdata(4,:) = sweepdata(4,:) + ch2r(:, ii)';
                    end % input channel 2 is microphone monitor of STIM...
                    Sweep_nn = Sweep_nn + 1;
                otherwise
            end
        else % not in alternation mode
            sweepdata(1,:) = sweepdata(1,:) +  ch1r(:, ii)';
            if(STIM.Monitor)
                sweepdata(3,:) = sweepdata(3,:) + ch2r(:, ii)';
            end
            Sweep_nonalt = Sweep_nonalt + 1;
        end
    end
    %fprintf(1, 'Ch2 (a) max in sweepdata: %f, nonalt: %d, last samp: %f\n', ...
    %    max(sweepdata(3,:)), Sweep_nonalt, max(ch2r(:, ii))');
    
    Total_np = Total_np + Sweep_np;
    Total_nn = Total_nn + Sweep_nn;
    Total_all = Total_all + Sweep_nonalt;
    %fprintf(1, 'i=%d, total_all = %d\n', i, Total_all);
    set(GUI.hrep, 'string', sprintf('%d', i*STIM.StimPerSweep)); % update acquisition trial counter
    if(mod(replot, 10) == 0)
        replot = 0;
    else
        replot = replot + 1;
    end
    % average the data from accepted sweeps, then sum in to the data arrays
    if(STIM.Alternate)
        sweepdata(1,:) = sweepdata(1,:)/Sweep_np;
        sweepdata(1,:) = sweepdata(1,:) - mean(sweepdata(1,tbs:tbe));
        data(1,:) = data(1,:) + sweepdata(1,:);
        sweepdata(2,:) = sweepdata(2,:)/Sweep_nn;
        sweepdata(2,:) = sweepdata(2,:) - mean(sweepdata(2,tbs:tbe));
        data(2,:) = data(2,:) + sweepdata(2,:);
        if STIM.Monitor
            sweepdata(3,:) = sweepdata(3,:)/Sweep_np;
            sweepdata(4,:) = sweepdata(4,:)/Sweep_nn;
            data(3,:) = data(3,:) + sweepdata(3,:);
            data(4,:) = data(4,:) + sweepdata(4,:);
        end
        set(PLOTS.data, 'clipping','On');
        hold on; % set(PLOT.PLOTS.data, 'hold', 'On');
        runningAvg = (data(1,:)/i+data(2,:)/i)/2.0;
        plotSweep{1,i} = plot(PLOTS.data, timebase_Display, (sweepdata(1,:)+sweepdata(2,:))/2.0, 'b-');
        plotSweep{1,i}.Color=[0.6, 0.3, 0.3, 0.5];
        
        if i == 1
            avgplot = plot(PLOTS.data, timebase_Display, runningAvg, 'k-', 'LineWidth', 2.0);
        else
            set(avgplot,'Xdata',timebase_Display,'YData',runningAvg);
        end
        % plot(PLOTS.data, timebase_Display, sweepdata(1,:), 'g-');
        % plot(PLOTS.data, timebase_Display, sweepdata(2,:), 'r-');
        md = max(abs(sweepdata(1,:)));
        
    else
        sweepdata(1,:) = sweepdata(1,:)/Sweep_nonalt;
        sweepdata(1,:) = sweepdata(1,:) - mean(sweepdata(1,tbs:tbe));
        data(1,:) = data(1,:) + sweepdata(1,:);
        runningAvg = data(1,:)/i;
        if STIM.Monitor
            sweepdata(3,:) = sweepdata(3,:)/Sweep_nonalt;
            data(3,:) = data(3,:) + sweepdata(3,:);
        end
        
        plotSweep{1, i} = plot(PLOTS.data, timebase_Display, sweepdata(1,:), 'b-');
        if i == 1
            avgplot = plot(PLOTS.data, timebase_Display, runningAvg, 'k-', 'LineWidth', 2.0);
        else
            set(avgplot,'Xdata',timebase_Display,'YData',runningAvg);
        end
        plotSweep{1,i}.Color=[0.6, 0.6, 0.6, 0.5];
        %         uistack(plotSweep{1,i}, 'top');
        md = max(abs(sweepdata(1,:)));
    end
    % trace_scale(md, PLOTS.data);
    % set(PLOTS.data, 'YLim', [-5e-6, 5e-6]);
    set_scale(PLOTS)
    %     ds = findobj('Tag', 'ABR_DisplayScale');
    %     dsstr = get(ds, 'String');
    %     dsval = get(ds, 'Value');
    %     if strcmp(dsstr{dsval}, 'auto') == 0
    %         dsscale = str2double(dsstr{dsval});
    %         set(PLOTS.data, 'YLim', [-dsscale dsscale]);
    %     end;
    set(PLOTS.data, 'xlim', [0 STIM.avg_dur]);
    set(PLOTS.data, 'clipping', 'on');
    drawnow;
    
    % fprintf(1, 'Ch2 max in sweepdata: %f, nonalt: %d\n', max(sweepdata(3,:)), Sweep_nonalt);
    if(STIM.Monitor == 1) % only do this if we have the monitor set on
        hold off
         fprintf(1, 'plotting max: %f  x %f\n', max(sweepdata(3,:)), max(timebase_Display));
        plot(PLOTS.signal2, timebase_Display, sweepdata(3,:), 'color', 'red');
        if(STIM.Alternate)
            hold on
            plot(PLOTS.signal2, timebase_Display, sweepdata(4,:), 'color', 'cyan');
        end
        md = max(abs(sweepdata(3,:)));
        set(PLOTS.signal2, 'xLimMode', 'Auto');
        set(PLOTS.signal2, 'yLimMode', 'Auto');
        %trace_scale(md, PLOTS.signal2);
    end
    
    % delay to next sweep:
    if i < STIM.NSweeps % (but not on last sweep)
        while toc < STIM.InterSweepInterval
            [HW, stf] = check_stop(HW, 0);
            if(stf == 1) % while waiting, check for stop.
                fprintf(2, 'Sweep interval stop\n');
                return;
            end
        end
    end
end  % END of the big FOR loop

% divde accumulated averages by the number of sweeps presented.

data(1,:) = data(1,:)/STIM.NSweeps;
data(2,:) = data(2,:)/STIM.NSweeps;
if STIM.Monitor
    data(3,:) = data(3,:)/STIM.NSweeps;
    data(4,:) = data(4,:)/STIM.NSweeps;
end

% figure(10);
% plot(data(1,:));
STIM.ACQPars_np = Total_np;
STIM.ACQPars_nn = Total_nn;
STIM.ACQPars_nall = Total_all;
STIM.ACQPars_tb = timebase_Display;
set_status('Done');
% fid = fopen('Tessaraw.txt','w');
% fprintf(fid,'%6.2f  %12.8f\n',rawdata);
% fprintf(2, 'max ch 2: %f\n', md);

end

