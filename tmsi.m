function score = tmsi(Xnew,Xtemplate,tau,r)
% Xnew [samples, channels]
% Xtemplate [freq, samples, channels]
% tau: temporally local range, constant
% r: weight for Tukeys tricube weighting function, constant
% Origin paper setting: tau = 24, r = 3;
if nargin == 2
    tau = 24; r = 3;
end
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
    Y = squeeze(Xtemplate(freq,:,:));
    Z = [Xnew';Y']';
    W = TukeysWeight(Nsample,tau,r);
    D = diag(sum(W,2));
    L = (D - W)';
    C = 1/Nsample*(Z'*L*Z);
    Cx = C(1:dimOne,1:dimOne);
    Cy = C(dimOne+1:dimA,dimOne+1:dimA);
    U = [inv(sqrtm(Cx)) zeros(dimOne,dimTwo);...
        zeros(dimTwo,dimOne) inv(sqrtm(Cy))];
    R = U*C*U';
    lamda = eig(R);
    nlamda = lamda/sum(lamda);
    score(freq) = 1 + sum(bsxfun(@times,nlamda,log(nlamda)))/log(dimA);
end
end

function W = TukeysWeight(Nsample,tau,r)
W = zeros(Nsample);
for i = 1:Nsample
    for j = 1:Nsample
        v = (j - i)/tau;
        if abs(v) < 1
            W(i,j) = (1 - abs(v)^r)^r;
        else
            W(i,j) = 0;
        end
    end
end
end
