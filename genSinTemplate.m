function template = genSinTemplate(stimuFreq, fsample, time, Nharmonic)
template = zeros(length(stimuFreq),fsample*time,2*Nharmonic);
timeSeq = 0:1/fsample:time-1/fsample;
for freqIdx = 1:length(stimuFreq)
    for i = 1:Nharmonic
        template(freqIdx,:,(i-1)*2+1) = sin(2*pi*i*stimuFreq(freqIdx)*timeSeq);
        template(freqIdx,:,i*2) = cos(2*pi*i*stimuFreq(freqIdx)*timeSeq);
    end
end

% Generate artifical sinosudial template for conventional CCA
% stimuFreq: Frequencies we use in visual stimulus
% fsample: Sampling frequency
% time: Time length of the signal (= Time window length)
% Nharmonic: Number of harmonic we consider in the template