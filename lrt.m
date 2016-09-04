function score = lrt(Xnew,Xtemplate)
% Xnew [samples,channels]
% Xtemplate [freq,samples,channels]
if ismatrix(Xtemplate)
    Xtemplate = reshape(Xtemplate,[1 size(Xtemplate)]);
end
score = zeros(size(Xtemplate,1),1);
Nsample = size(Xnew,1);
Xnew = normalizeSignal(Xnew,1);
Xtemplate = normalizeSignal(Xtemplate,2);
dimOne = size(Xnew,2);
dimTwo = size(Xtemplate,3);
dimA = dimOne + dimTwo;
for freq = 1:size(Xtemplate,1)
    subRef = squeeze(Xtemplate(freq,:,:));
    intSignal = [Xnew';subRef']';
    C = 1/Nsample*(intSignal'*intSignal);
    Cx = C(1:dimOne,1:dimOne);
    Cy = C(dimOne+1:dimA,dimOne+1:dimA);
    score(freq) = 1 - (det(C)/(det(Cx)*det(Cy)+10e-9))^(1/size(Xtemplate,3));
end
end

% Reference:
% Zhang, Y., Dong, L., Zhang, R., Yao, D., Zhang, Y., Xu, P.: An Efficient 
% Frequency Recognition Method Based on Likelihood Ratio Test for SSVEP-Based 
% BCI. Computational and Mathematical Methods in Medicine. 2014, 1¡V7 (2014). 
