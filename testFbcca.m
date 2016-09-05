clear;
clc;
addpath('C:\SSVEP\Algorithm')
%%  Data preprocessing: Janir Database
filename = 'ReadMe (dataset_processed) english.txt';
[stimuFreq,fsample,dataName] = extrectReadme(filename);
stop
fprintf('Stimulus frequency: \n');
disp(stimuFreq);
freqLength = length(stimuFreq);
time = 4; %second
numSubject = length(dataName);
Nharmonic = 3;
ridx = 1;

% Data preprocessing (Read -> Normalize -> reshape -> concatenate)
for i = 1:numSubject
    temp = load(dataName{i});
    data = temp.ssvep_data;
    data = data(:,:,:,1:time*fsample);
    if ~exist('allData','var')
        sampleLength = size(data,4);
        channelLength = size(data,3);
        trialLength = size(data,2);
        allData = cell(numSubject,1);
    end
    data = permute(data,[2 1 4 3]);
    allData{i} = data;
end
%%  Generate  template
% Generate sinosudial template
sinTemplate = genSinTemplate(stimuFreq, fsample, time, Nharmonic);

% Generate average template for every subject's data
sourceTemplate = cell(numSubject,1);
for subject = 1:numSubject
    ssvep = allData{subject};
    sourceTemplate{subject} = squeeze(mean(ssvep,1));
end
%%
rec = zeros(numSubject,1);
for subject = 1:numSubject
    ssvep = allData{subject};
    for trial = 1:trialLength
        for freq = 1:freqLength
            X = squeeze(ssvep(trial,freq,:,:));
            p =  mec(X, sinTemplate, 0.1, fsample, Nharmonic, 0.25);
            [~,maxLoc] = max(p);
            if maxLoc == freq, rec(subject) = rec(subject) + 1; end;
        end
    end
end
acc = rec/(trialLength*freqLength)*100
stop

% startFreq = min(stimuFreq);
% endFreq = max(stimuFreq);
% stopFreq = 60;
% Nfilter = 7;
% tol = 2;
% filterBank = genFB(startFreq,endFreq,stopFreq,fsample,Nfilter,tol);
% 
% ss = 0:0.25:2;
% for aa = 3:9
%     for bb = 1:5
% a = ss(aa);
% b = ss(bb);
% rec = zeros(numSubject,1);
% w = [1:Nfilter].^-a + b;
% for subject = 1:numSubject
%     ssvep = allData{subject};
%     for trial = 1:trialLength
%         for freq = 1:freqLength
%             X = squeeze(ssvep(trial,freq,:,:));
%             p = zeros(freqLength,Nfilter);
%             for sb = 1:Nfilter
%                 Xsb = filtfilt(filterBank{sb},X);
%                 p(:,sb) = ccaExtend(Xsb,[], sinTemplate, 'CCA');
%             end
%             fp = bsxfun(@times,w,p.^2);
%             [~,maxLoc] = max(sum(fp,2));
%             if maxLoc == freq, rec(subject) = rec(subject) + 1;end; 
%         end
%     end
% end
% 
% fprintf('a = %d, b = %d: \n',a,b)
% acc = rec/(trialLength*freqLength);
% disp(mean(acc))
%     end
% end