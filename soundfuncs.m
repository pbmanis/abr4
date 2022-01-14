classdef soundfuncs
%{
This class provides functions for dealing with sound levels and 
octave calculations
compute_spl, dBPerVPa and mVPerPa compute the sound pressure level or
calibration factors from a microphone voltage.
    MIC is a structure with the following fields:

          Gain: 20  # amplifier gain
        RefSig: 94  # level used for the calibration in this file
          Vrms: 0.0395  # raw rms microphone voltage, in V
        Vref_c: 0.0286  # cosinor measure of microphone voltage, in Vrms
       Vref_bp: 0.0389  # bandpassed measure of microphone voltae, in Vrms
    Microphone: '7016#10252'  # identiy of the the microphone
          Date: '11-Jan-2022'  # date of calibration
      dBPerVPa: -48.2029   # sensitivity
       mVPerPa: 3.8892   # calibration factor

octave_calc computes the 1/n'th octave frequencies above and below a
    center frequency.


%}


    methods(Static)
        function [dBSPL] = compute_spl(mic_vrms_V, MIC)
            % compute_spl from a microphone voltage (rms, in Volts), given the
            % MIC structure that has calibration information in it (mic amplifier gain,
            % reference calibration dB, and microphone sensitivity.
            dBSPL = MIC.RefSig - MIC.Gain + 20.0*log10(mic_vrms_V*1e3/MIC.mVPerPa);
            
        end
        
        function [dBPerVPa] = compute_dBPerVPa(mic_vrms_V, MIC)
            % compute the microphone sensitivity given the microphone
            % output rms voltage at a given frequency in Volts.
            % sensitivity is usually stated relative to 94 dB SPL
            dBPerVPa = 20.0*log10(mic_vrms_V) - MIC.Gain - (MIC.RefSig - 94.0); % sensitivity stated relative to 1Pa.
            
        end
        function [mVPerPa] = compute_mVPerPa(MIC)
            mVPerPa = 1000.0*10.0^(MIC.dBPerVPa/20.0);
        end
        
        function [bp_freqs] = octave_calc(f0, oct, Fs)
            % compute the low and high frequencies corresponding to 
            % the octave fraction. Limit the values according to the
            % nyquist limit based on the sampling ferquency Fs.
            % for example, with oct = 3, the two frequences form the 
            % 1/3 octave filter frequencies.
            octfrac = power(2.0, 1.0/oct);  % fraction of an octave filter
            flow = f0/octfrac;
            fhigh = f0*octfrac;
            if fhigh > Fs/2.0
                fhigh = floor(Fs/2.0);
            end
            bp_freqs = [flow, fhigh];
        end
        
        function [spl] = spl_at_f(freqs, spls, f)
            % compute the spl at frequencies in f by doing a linear
            % interpolation from the measured freq vs spl plot
            % from the calibration
            spl = interp1(log2(freqs), spls, log2(f), 'linear');
        end
    end
end


