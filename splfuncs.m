%{
This class provides functions that compute the sound pressure level or
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

%}

classdef splfuncs
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
            octfrac = power(2.0, 1.0/oct);  % fraction of an octave filter
            flow = f0/octfrac;
            fhigh = f0*octfrac;
            if fhigh > Fs/2.0
                fhigh = floor(Fs/2.0);
            end
            bp_freqs = [flow, fhigh];
        end
    end
end


