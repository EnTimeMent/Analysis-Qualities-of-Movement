function [Sf] = LowPassButtDouble (S, SampFreq, CutFreq)
    % LowPassButtDouble : low pass filter of a signal with dual pass Butterworh
    %
    % Calling Sequence
    %  [Sf] = LowPassButtDouble (S, SampFreq, CutFreq)
    %
    % Parameters
    %  S        : vector,  the input signal S must be vector.
    %  SampFreq : number,  sampling frequency in Hz
    %  CutFreq  : number,  cutoff frequency in Hz
    %  Sf       : vector,  the signal after the filtering
    %
    % Description
    %  LowPassButtDouble : low pass filter of a signal with dual pass Butterworh
    %  to cancel phase shift. The resulting filter has an order of 4 and an
    %  effective cut frequency of about CutFreq * 83%.
    %  Sf is a vector (same size as S)
    %
    % Examples
    %  To get a signal filtered at 8Hz
    %      T = linspace(0,10,1001);
    %      S = sin(2*%pi*T)+sin(20*%pi*T);
    %      Sf = LowPassButtDouble (S, SampFreq, 8);
    %      plot(T, S, '-k', T, Sf, '-b')
    %
    % Authors
    %  Denis Mottet - Univ Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- D. Mottet -- 2016-05-28
    %    First version
    %  Version 1.0.1 -- D. Mottet -- 2018-03-27
    %    Documentation following
    %      https:%wiki.scilab.org/Guidelines%20To%20Design%20a%20Module
    %  Version 1.0.2 -- D. Mottet -- 2019-10-09
    %      self contained example (cut-paste is ok)



    if min(size(S)) > 1 then
        error('The input signal should be a vector')
    end

    % Correction of cut frequency du to the dual pass in fltsflts
    Cutoff = CutFreq ./ 0.802 ;
    % This is a simplification (TODO : improve this function)
    % https:%www.codeproject.com/Articles/1267916/Multi-pass-Filter-Cutoff-Correction
    
    % Computation of the Butterworth filter equation
    [b,a] = butter(2,Cutoff./SampFreq);
    
    % Dual-pass filtering in the time domain
    Sf = filtfilt(b, a, S);

end
