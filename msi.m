function p = msi(Xnew,Xtemplate)
% Xnew [samples, channels]
% Xtemplate [freq, samples, channels]
if ismatrix(Xtemplate)
    Xtemplate = reshape(Xtemplate,[1 size(Xtemplate)]);
end
p = zeros(size(Xtemplate,1),1);
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
    U = [inv(sqrtm(Cx)) zeros(dimOne,dimTwo);...
        zeros(dimTwo,dimOne) inv(sqrtm(Cy))];
    R = U*C*U';
    lamda = eig(R);
    nlamda = lamda/sum(lamda);
    p(freq) = 1 + sum(bsxfun(@times,nlamda,log(nlamda)))/log(dimA);
end
end

% Zhang, Y., Xu, P., Cheng, K., Yao, D.: Multivariate synchronization index 
% for frequency recognition of SSVEP-based brain¡Vcomputer interface. Journal 
% of Neuroscience Methods. 221, 32¡V40 (2014). 
