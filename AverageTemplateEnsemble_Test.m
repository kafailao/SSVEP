clear;clc;
%% Testing parameter
time = 1;
Nh = 1;
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

% Calculate the average EEG response template for each subject
allTemplate = cell(numSubject,1);
for subject = 1:numSubject
    allTemplate{subject} = squeeze(mean(allData{subject}));
end

% One-hot encoded label for each data
label = zeros(numSubject - 1,trialLength*freqLength,freqLength);
for i = 1:numSubject
    label(i,:,:) = repmat(eye(freqLength),trialLength,1);
end

% Artifical(sinusoidal template)
sinTemplate = genSinTemplate(stimuFreq,fsample,time,Nh);
%% Extract features for all data (including test data)
startIdx = round(fsample*0.135); %Cut the first 0.135 second
nF = 2; % Number of features estimate from each template
featureSet = zeros(numSubject,trialLength*freqLength,numSubject*nF+1,freqLength);
% Train generic ensemble classifier
for testSubject = 1:numSubject
    ssvep = allData{testSubject};
    for templateSubject = 1:numSubject
        Xtemplate = allTemplate{templateSubject};
        Xtemplate = Xtemplate(:,startIdx+1:startIdx+time*fsample,:);
        for trial = 1:trialLength
            for freq = 1:freqLength
                Xnew = squeeze(ssvep(trial,freq,startIdx+1:startIdx+time*fsample,:));
                [~,pvec] = ccaExtend(Xnew,Xtemplate,sinTemplate,'Combination3');
                if nF == 1
                    p = sign(pvec(:,1)).*pvec(:,1).^2 + sign(pvec(:,3)).*pvec(:,3).^2;
                    featureSet(testSubject,(trial - 1)*freqLength + freq,templateSubject+1,:) = p;
                elseif nF == 2
                    featureSet(testSubject,(trial - 1)*freqLength + freq,(templateSubject-1)*nF+2,:) = sign(pvec(:,1)).*pvec(:,1).^2;
                    featureSet(testSubject,(trial - 1)*freqLength + freq,templateSubject*nF+1,:) = sign(pvec(:,3)).*pvec(:,3).^2;
                end
                if templateSubject == 1
                    featureSet(testSubject,(trial - 1)*freqLength + freq,1,:) = sign(pvec(:,2)).*pvec(:,2).^2;
                end
            end
        end
    end
end
ppath = 'C:\SSVEP\Data\';
str = sprintf('Janir_Nh%d_minMSE_%ds',Nh,time);
if exist(ppath,'file') == 0, mkdir(ppath); end;
save([ppath str '.mat'],'featureSet','numSubject','subjectSeq','nF','label','trialLength','freqLength');

rec = zeros(numSubject,1);
lassoApproach = 'IndexMinMSE';
Beta = zeros((numSubject-1)*nF+1,freqLength);
for targetSubject = 1:numSubject
    %% Use LASSO to combine feature toward 
    sourceSubject = find(subjectSeq ~= targetSubject);
    trainSet = featureSet(sourceSubject,:,:,:);
    if nF == 1
        trainSet(:,:,targetSubject+1,:) = [];
    else
        trainSet(:,:,(targetSubject-1)*nF+2:targetSubject*nF+1,:) = [];
    end
    for freq = 1:freqLength
        rTrainSet = reshape(squeeze(trainSet(:,:,:,freq)),size(trainSet,1)*size(trainSet,2),size(trainSet,3));
        rLabel = reshape(squeeze(label(sourceSubject,:,freq)),(size(label,1)-1)*size(label,2),1);
        [B,FitInfo] = lasso(rTrainSet,rLabel,'CV',10);
        Beta(:,freq) = B(:, FitInfo.IndexMinMSE);
%         Beta(:,freq) = B(:,FitInfo.Index1SE);
    end
    
    %% Recognition
    coef = squeeze(featureSet(targetSubject,:,:,:));
    if nF == 1
        coef(:,targetSubject+1,:) = [];
    else
        coef(:,(targetSubject-1)*nF+2:targetSubject*nF+1,:) = [];
    end
    
    fcoef = zeros(freqLength,1);
    for trial = 1:trialLength
        for freq = 1:freqLength
            for tFreq = 1:freqLength
                fcoef(tFreq) = squeeze(coef((trial-1)*freqLength + freq,:,tFreq))*Beta(:,tFreq);
            end
            [~,maxLoc] = max(fcoef);
            if maxLoc == freq, rec(targetSubject) = rec(targetSubject) + 1;end
        end
    end
    fprintf('Accuracy of S%d = %.2f\n',targetSubject,rec(targetSubject)/(trialLength*freqLength));
end
acc = zeros(numSubject+2, 1);
acc(1:numSubject) = rec/(trialLength*freqLength)*100;
acc(numSubject+1) = mean(acc(1:numSubject));
acc(numSubject+2) = std(acc(1:numSubject));
%%  Save result to xls
% col_header = strsplit(num2str(timeSeq));     %Row cell array (for column labels)
col_header = {'1s'};
row_header = cell(numSubject + 2,1);
for i = 1:numSubject
    row_header{i}=['s' num2str(i)];     %Column cell array (for row labels)
end
row_header{numSubject+1} = 'Mean';
row_header{numSubject+2} = 'Std';
xlswrite('ATE_JanirDataset.xlsx',acc,sprintf('Nh = %d',Nh),'B2');     %Write data
xlswrite('ATE_JanirDataset.xlsx',col_header,sprintf('Nh = %d',Nh),'B1');     %Write column header
xlswrite('ATE_JanirDataset.xlsx',row_header,sprintf('Nh = %d',Nh),'A2');      %Write row header