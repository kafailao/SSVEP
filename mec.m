function p = mec(Xnew, Xtemplate, lamda, fsample, Nharmonic, alpha)
% [Input]
% Xnew [samples,channels]: Test EEG data
% Xtemplate [freq,samples,channels]: SSVEP related components, can be
% sinosudial template or EEG template
% lamda: Ratio of noise energy retain
% fsample: Sampling frequency
% Nharmonic: Number of harmonic contains in Xtemplate, use 1 if EEG
% template is used
% alpha: Paramter used in calculating softmax probaility, alpha = 0.25
time = 0:1/fsample:size(Xnew,1)/fsample - 1/fsample;
P = zeros(size(Xtemplate,1),1);
% Remove power line interference
Xpower = [sin(2*pi*50*time);cos(2*pi*50*time);...
    sin(2*pi*60*time);cos(2*pi*60*time)]';
Xnew = Xnew - Xpower*inv(Xpower'*Xpower)*Xpower'*Xnew;

for freq = 1:size(Xtemplate,1)
    Xssvep = squeeze(Xtemplate(freq,:,:));
    % Remove potential SSVEP components at ith frequency
    Xnoise = Xnew - Xssvep*inv(Xssvep'*Xssvep)*Xssvep'*Xnew;
    [eigVector,~,eigValue] = princomp(Xnoise);
    idx = cumsum(eigValue)/sum(eigValue) > (1 - lamda);
    eigVector = bsxfun(@rdivide,eigVector,sqrt(eigValue)');
    W = eigVector(:,idx);
    s = Xnew*W;
    recP = 0;
    for channel = 1:size(s,2)
        for Nh = 1:Nharmonic
            recP = recP + norm(Xssvep(:,2*(Nh - 1)+1:2*Nh)'*s(:,channel),2)^2;
        end
    end
    P(freq) = recP;
end

pNorm = P/sum(P);
pSoftmax = exp(alpha*pNorm)/sum(exp(alpha*pNorm));
p = pSoftmax;

% [Orginal from]
% Friman, O., Volosyak, I., Graser, A.: Multiple Channel Detection of Steady-
% State Visual Evoked Potentials for Brain-Computer Interfaces. IEEE Transac-
% tions on Biomedical Engineering IEEE Trans. Biomed. Eng. 54, 742¡V750 (2007). 

% [Modified as]
% Volosyak, I.: SSVEP-based Bremen¡VBCI interface¡Xboosting information 
% transfer rates. J. Neural Eng. Journal of Neural Engineering. 8, 036020 (2011). 
