function Xnoise = getNoisePattern(Xnew,signalFreq, fsample, Nharmonic)
% [Input]
% Xnew [samples,channels]: Test EEG data
% signalFreq: stimuli frequency that evoked Xnew
% fsample: Sampling frequency, constant
% Nharmonic: Number of harmonic contains in Xtemplate, use 1 if EEG

if isscalar(signalFreq)
    Xssvep = squeeze(genSinTemplate(signalFreq,fsample,size(Xnew,1)/fsample,Nharmonic));
elseif ismatrix(signalFreq)
    Xssvep = signalFreq;
end
    
time = 0:1/fsample:size(Xnew,1)/fsample - 1/fsample;
% Remove power line interference
Xpower = [sin(2*pi*50*time);cos(2*pi*50*time);...
    sin(2*pi*60*time);cos(2*pi*60*time)]';
Xnew = Xnew - Xpower*inv(Xpower'*Xpower)*Xpower'*Xnew;

% Remove potential SSVEP components at ith frequency
Xnoise = Xnew - Xssvep*inv(Xssvep'*Xssvep)*Xssvep'*Xnew;