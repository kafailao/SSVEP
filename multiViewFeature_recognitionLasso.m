clear;clc;
NhSeq = 1:3;
timeSeq = 0.5:0.5:3.5;
nameDataset = 'JanirDataset';
method = 'lasso';
subApproach = 'IndexMinMSE';
numSubject = 28;
subjectSeq = 1:numSubject;
freqLength = 10;
trialLength = 5;
nF = 2;

rec_BetaRatio = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength);
rec_confusion = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength,freqLength);
for Nhidx = 1:length(NhSeq)
    Nh = NhSeq(Nhidx);
    acc = zeros(numSubject+2,length(timeSeq));
    for tidx = 1:length(timeSeq)
        time = timeSeq(tidx);
        load(['Feature\' sprintf('ECCA3_%s_Nh%d_time%d.mat',nameDataset,Nh,time*10)]);        
        rec = zeros(numSubject,1);
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
            rec_BetaRatio(Nhidx,tidx,targetSubject,:) = Beta(1,:)/sum(abs(Beta));
            
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
                    rec_confusion(Nhidx,tidx,targetSubject,freq,maxLoc) = rec_confusion(Nhidx,tidx,targetSubject,freq,maxLoc) + 1;
                end
            end
            fprintf('Accuracy of S%d = %.2f\n',targetSubject,rec(targetSubject)/(trialLength*freqLength));
        end
        acc(1:numSubject,tidx) = rec/(trialLength*freqLength)*100;
        acc(numSubject+1,tidx) = mean(acc(1:numSubject,tidx));
        acc(numSubject+2,tidx) = std(acc(1:numSubject,tidx));
    end
    %%  Save result to xls
    col_header = strsplit(num2str(timeSeq));     %Row cell array (for column labels)
    col_header = strcat(col_header,'s');
%     col_header = {sprintf('%1.1f',time)};
    row_header = cell(numSubject + 2,1);
    for i = 1:numSubject
        row_header{i}=['s' num2str(i)];     %Column cell array (for row labels)
    end
    row_header{numSubject+1} = 'Mean';
    row_header{numSubject+2} = 'Std';
    filename = ['Result\' method '_' subApproach '_' nameDataset '.xlsx'];
    xlswrite(filename,acc,sprintf('Nh = %d',Nh),'B2');     %Write data
    xlswrite(filename,col_header,sprintf('Nh = %d',Nh),'B1');     %Write column header
    xlswrite(filename,row_header,sprintf('Nh = %d',Nh),'A2');      %Write row header
end
save(['Result\' method '_' subApproach '_' nameDataset '_confusion.mat'],'rec_confusion','rec_BetaRatio');