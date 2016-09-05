clear;clc;
%% Testing parameter
% nameDataset = 'UCSDDataset';
nameDataset = 'JanirDataset';
filterOn = false;
%% Initialization
% addpath('C:\SSVEP\Algorithm')
loadDataTime = 4;
readmeFileName = 'readme.txt';

% Load data to workspace
% allData: SSVEP data in cell format
% allData{j} contains SSVEP data from j^th subject
[allData,stimuFreq,fsample,dataSize] = prepareData(loadDataTime,nameDataset,readmeFileName,filterOn);

% Basic info of data
trialLength = dataSize(1); %Number of recorded EEG response for each stimulus frequency
freqLength = dataSize(2); %Number of visual stimulus 
sampleLength = dataSize(3); %Number of time points in each record
channelLength = dataSize(4); %Number of channels used in each experiment
numSubject = length(allData); %Number of subjects in the data set
trialSeq = 1:trialLength;
subjectSeq = 1:numSubject;

% Generate sinusoidal template
sinTemplate = genSinTemplate(stimuFreq, fsample, time, Nharmonic);
%%
startFreq = min(stimuFreq);
endFreq = max(stimuFreq);
stopFreq = 60;
Nfilter = 7;
tol = 2;
filterBank = genFB(startFreq,endFreq,stopFreq,fsample,Nfilter,tol);

ss = 0:0.25:2;
for aa = 3:9
    for bb = 1:5
a = ss(aa);
b = ss(bb);
rec = zeros(numSubject,1);
w = [1:Nfilter].^-a + b;
for subject = 1:numSubject
    ssvep = allData{subject};
    for trial = 1:trialLength
        for freq = 1:freqLength
            X = squeeze(ssvep(trial,freq,:,:));
            p = zeros(freqLength,Nfilter);
            for sb = 1:Nfilter
                Xsb = filtfilt(filterBank{sb},X);
                p(:,sb) = ccaExtend(Xsb,[], sinTemplate, 'CCA');
            end
            fp = bsxfun(@times,w,p.^2);
            [~,maxLoc] = max(sum(fp,2));
            if maxLoc == freq, rec(subject) = rec(subject) + 1;end; 
        end
    end
end

fprintf('a = %d, b = %d: \n',a,b)
acc = rec/(trialLength*freqLength);
disp(mean(acc))
    end
end