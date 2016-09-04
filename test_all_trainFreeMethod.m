clearvars -except method;clc;
%% Testing parameter
% time = 1;
% Nh = 1;
% method = 'MEC_AR';
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

%% Recognition
timeSeq = 0.5:0.5:3.5;
NhSeq = 1:3;
for Nhidx = 1:length(NhSeq)
    Nh = NhSeq(Nhidx);
    acc = zeros(numSubject+2,length(timeSeq));
    for tidx = 1:length(timeSeq)
        time = timeSeq(tidx);   
        % Artifical(sinusoidal template)
        sinTemplate = genSinTemplate(stimuFreq,fsample,time,Nh);
        startIdx = 0; %%%%%%%%%%%%%%%%%%5
        
        rec = zeros(numSubject,1);
        for targetSubject = 1:numSubject
            ssvep = allData{targetSubject};
            for trial = 1:trialLength
                for freq = 1:freqLength
                    Xnew = squeeze(ssvep(trial,freq,startIdx+1:startIdx+time*fsample,:));
                    switch method
                        case 'MEC_AR'
                            lambda = 0.1; p = 10; alpha = 0;
                            score = mec(Xnew, sinTemplate, fsample, Nh, lambda, stimuFreq, p, alpha);
                        case 'MEC1'
                            lambda = 0.1; p = 0; alpha = 0.25;
                            score = mec(Xnew, sinTemplate, fsample, Nh, lambda, stimuFreq, p, alpha);                            
                        case 'CCA'
                            score = ccaExtend(Xnew,sinTemplate,sinTemplate,'CCA');
                        case 'MCC_AR'
                            p = 10; alpha = 0;
                            score = mcc(Xnew, sinTemplate, fsample, Nh, stimuFreq, p, alpha);
                        case 'MCC1'
                            p = 0; alpha = 0.25;
                            score = mcc(Xnew, sinTemplate, fsample, Nh, stimuFreq, p, alpha);                            
                        case 'CVARS'
                            p = 10; alpha = 0;
                            score = cvars(Xnew, sinTemplate, fsample, Nh, stimuFreq, p, alpha);
                        case 'MSI'
                            score = msi(Xnew, sinTemplate);
                        case 'TMSI'
                            tau = 24; r = 3;
                            score = tmsi(Xnew,sinTemplate,tau,r);
                        case 'LASSO'
                            lambda = 0.1;
                            score = slasso(Xnew,sinTemplate,lambda);
                        case 'LRT'
                            score = lrt(Xnew, sinTemplate);
                        case 'DPLS'
                            score = dpls(Xnew, sinTemplate);
                    end     
                    [~,maxLoc] = max(score);
                    if maxLoc == freq, rec(targetSubject) = rec(targetSubject) + 1; end;
                end
            end
            fprintf('Accuracy of S%d: %.2f\n',targetSubject,rec(targetSubject)/(trialLength*freqLength));
        end
        acc(1:numSubject,tidx) = rec/(trialLength*freqLength)*100;
        acc(numSubject+1,tidx) = mean(acc(1:numSubject,tidx));
        acc(numSubject+2,tidx) = std(acc(1:numSubject,tidx));
    end
    %%  Save result to xls
    col_header = strsplit(num2str(timeSeq));     %Row cell array (for column labels)
    col_header = strcat(col_header,'s');
    row_header = cell(numSubject + 2,1);
    for i = 1:numSubject
        row_header{i}=['s' num2str(i)];     %Column cell array (for row labels)
    end
    row_header{numSubject+1} = 'Mean';
    row_header{numSubject+2} = 'Std';
    filename = ['Result\' method '_' nameDataset '.xlsx'];
    xlswrite(filename,acc,sprintf('Nh = %d',Nh),'B2');     %Write data
    xlswrite(filename,col_header,sprintf('Nh = %d',Nh),'B1');     %Write column header
    xlswrite(filename,row_header,sprintf('Nh = %d',Nh),'A2');      %Write row header
end