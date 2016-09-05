function filterBank = genFB(startFreq,endFreq,stopFreq,fsample,Nfilter,tol)
% startFreq: minimum pass frequency (e.g., the minimum stimulus frequency)
% endFreq: end frequency (e.g., the maximum stimulus frequency)
% stopFreq: maximum frequency for the filter bank
% fsample: sampling frequency
% Nfilter: number of filter bank
% tol: reduce/increase startFreq/stopFreq with tol

startFreq = floor(startFreq);
bandwidth = endFreq - startFreq;
filterBank = cell(Nfilter,1);
fs = fsample/2;
for n = 1:Nfilter
    filterBank{n} = designfilt('bandpassiir','FilterOrder',10, ...
         'HalfPowerFrequency1',startFreq + (n-1)*bandwidth - tol,'HalfPowerFrequency2',stopFreq + tol, ...
         'SampleRate',fsample);
end