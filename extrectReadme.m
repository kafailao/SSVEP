function [stimuFreq,fsample,dataName] = extrectReadme(filename)
fid = fopen(filename,'r','n','Unicode');
str = fread(fid,'*char')';
tmpcell = regexp(str,'(?<=Stimulus frequency:[^0-9]*)[0-9.\s]+','match');
stimuFreq = str2num(tmpcell{1});
tmpcell = regexp(str,'(?<=Sampling frequency:[\s]*)[0-9]+','match');
fsample = str2num(tmpcell{1});
[~,endIndex] = regexp(str,'Data name:');
tmpcell = regexp(str,['(?<=^.{' num2str(endIndex) ',' num2str(length(str)) '})(?=[\s]*).+?(?=[\s])'],'match');
dataName = tmpcell(~(cellfun('length',tmpcell) <= 1));
end
