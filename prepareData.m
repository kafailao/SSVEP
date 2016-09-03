function [allData,stimuFreq,fsample,dataSize] = prepareData(time,nameDataset,readmeFileName,filterOn)
% [Input]
% time: Window length (in second) of EEG data
% [Output]
% allData: cell array that stores EEG data from subjects
%          e.g., EEG data from ith subject can be extracted by ssvep = allData{i}
%          ssvep [Number of trial,Number of frequency,Number of
%          samples,Number of channels]
% stimuFreq: Value of stimulus frequency [1,Number of frequency]
% fsample: Sampling frequency
%%  Data preprocessing: Janir's Dataset
if strcmp(nameDataset,'UM2014') || strcmp(nameDataset,'JanirDataset')
    addpath('C:\SSVEP\UM_Dataset_Janir_2014')
    % Read filename from readme document
    [stimuFreq,fsample,dataName] = extrectReadme(readmeFileName);
    fprintf('Stimulus frequency: \n');
    disp(stimuFreq);
    fprintf('Sampling frequency: \n');
    disp(fsample);    
    numSubject = length(dataName);
    
    % Filter design
    if filterOn
        bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
            'HalfPowerFrequency1',5,'HalfPowerFrequency2',60, ...
            'SampleRate',fsample);
    end
    
    % Data preprocessing (Read -> Normalize -> reshape -> concatenate)
    for i = 1:numSubject
        temp = load(strtrim(dataName{i}));
        data = temp.ssvep_data;
        data = data(:,:,:,1:time*fsample);
        if ~exist('allData','var')
            allData = cell(numSubject,1);
        end
        data = permute(data,[2 1 4 3]);
        if filterOn
            for trial = 1:size(data,1)
                for freq = 1:size(data,2)
                    data(trial,freq,:,:) = filtfilt(bpFilt,squeeze(data(trial,freq,:,:)));
                end
            end
        end
        allData{i} = data;
    end
    dataSize = size(data);
    rmpath('C:\SSVEP\UM_Dataset_Janir_2014');
end
%% Data preprocessing: UCSD's Dataset
% Nakanishi, M., Wang, Y., Wang, Y.-T., Jung, T.-P.: A Comparison Study of 
% Canonical Correlation Analysis Based Methods for Detecting Steady-State 
% Visual Evoked Potentials. PLOS ONE. 10, (2015). 
if strcmp(nameDataset,'UCSD') || strcmp(nameDataset,'UCSDDataset')
    addpath('C:\SSVEP\UCSD Dataset_2015');
    stimuFreq = [9.25, 11.25, 13.25, 9.75, 11.75, 13.75, 10.25, 12.25, 14.25, 10.75, 12.75, 14.75];
    fsample = 256;
    numSubject = 10;
    
    % Create data name 
    dataName = cell(numSubject,1);
    for subject = 1:numSubject
        dataName{subject} = ['s' num2str(subject) '.mat'];
    end
    fprintf('Stimulus frequency: \n');
    disp(stimuFreq);
    fprintf('Sampling frequency: \n');
    disp(fsample);  
    numSubject = length(dataName);
    
    % Filter design
    if filterOn
        bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
            'HalfPowerFrequency1',6,'HalfPowerFrequency2',80, ...
            'SampleRate',fsample);
    end
    
    % Extract SSVEP segment (Stimulus on set at 39th sample) 
    delay = 39 + round(0.135*fsample);
    
    % Data preprocessing (Read -> reshape -> concatenate)
    for i = 1:numSubject
        temp = load(strtrim(dataName{i}));
        data = temp.eeg;
        data = data(:,:,delay+1:delay+fsample*time,:);
        if ~exist('allData','var')
            allData = cell(numSubject,1);
        end
        data = permute(data,[4 1 3 2]);
        if filterOn
            for trial = 1:size(data,1)
                for freq = 1:size(data,2)
                    data(trial,freq,:,:) = filtfilt(bpFilt,squeeze(data(trial,freq,:,:)));
                end
            end
        end
        allData{i} = data;
    end
    dataSize = size(data);
    rmpath('C:\SSVEP\UCSD Dataset_2015');
end
end