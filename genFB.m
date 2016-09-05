function filterBank = genFB(startFreq,endFreq,stopFreq,fsample,Nfilter,tol)
startFreq = floor(startFreq);
bandwidth = endFreq - startFreq;
filterBank = cell(Nfilter,1);
fs = fsample/2;
for n = 1:Nfilter
%     fst = startFreq + (n-1)*bandwidth;
%     fed = stopFreq;
%     [N,Wn] = cheb1ord([fst fed]/fs,[fst - tol fed + tol]/fs,3,60);
%     [b,a] = cheby1(N,1,Wn);
%     filterBank{n,1} = b;
%     filterBank{n,2} = a;
    filterBank{n} = designfilt('bandpassiir','FilterOrder',10, ...
         'HalfPowerFrequency1',startFreq + (n-1)*bandwidth - tol,'HalfPowerFrequency2',stopFreq + tol, ...
         'SampleRate',fsample);
%     designfilt('bandpassiir','FilterOrder',10, ...
%         'PassbandFrequency1',startFreq + (n-1)*bandwidth - tol,'PassbandFrequency2',stopFreq + tol, ...
%         'PassbandRipple',1,'SampleRate',fsample);
end