clear;clc;
%% Testing parameter
nameDataset = 'JanirDataset';
filterOn = false;
nF = 2; % Number of features estimate from each template
%% Initialization
loadDataTime = 4;
readmeFileName = 'readme.txt';

% Load data to workspace
% allData: SSVEP data in cell format
% allData{j} contains SSVEP data from j^th subject
[allData,stimuFreq,fsample,dataSize] = prepareData(loadDataTime,nameDataset,readmeFileName,filterOn);
startIdx = round(fsample*0.135); %Cut the first 0.135 second, not applicable for 'UCSDDataset'

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

%% Loop through different time length and Nh
NhSeq = 1:3;
timeSeq = 0.5:0.5:3.5;

for Nhidx = 1:length(NhSeq)
    Nh = NhSeq(Nhidx);
    acc = zeros(numSubject+2,length(timeSeq));
    for tidx = 1:length(timeSeq);        
        time = timeSeq(tidx);
        % Artifical(sinusoidal template)
        sinTemplate = genSinTemplate(stimuFreq,fsample,time,Nh);
        %% Extract features for all data (including test data)
        % featureSet: [testSubject,size of data, views, features]
        % Each subject provide nF views on the input data + one view from
        % CCA with sinusodial template
        featureSet = zeros(numSubject,trialLength*freqLength,numSubject*nF+1,freqLength);
        % Extract feature through CCA with artifical/subject-specific
        % template
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
        save(['Feature\' sprintf('ECCA3_%s_Nh%d_time%d.mat',nameDataset,Nh,time*10)],'featureSet','label');
        fprintf(sprintf('Have saved ECCA3_%s_Nh%d_time%d.mat\n',nameDataset,Nh,time*10));
    end
end