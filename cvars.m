function score = cvars(Xnew, Xtemplate, fsample, Nharmonic, stimuFreq, p, alpha)
% [Input]
% Xnew [samples,channels]: Test EEG data
% Xtemplate [freq,samples,channels]: SSVEP related components, can be
% sinosudial template or EEG template
% fsample: Sampling frequency, constant
% Nharmonic: Number of harmonic contains in Xtemplate, use 1 if EEG
% template is used
% lamda: Ratio of noise energy retain
% stimuFreq: frequency used for stimulus
% p: Order of autoregressive model for estimate noise
% if p = 0, assume noise power = 1 (consider signal power only)
% alpha: Paramter used in calculating softmax probaility, alpha = 0.25
% Leave alpha = 0 if softmax is not used

time = 0:1/fsample:size(Xnew,1)/fsample - 1/fsample;
score = zeros(size(Xtemplate,1),1);
% Remove power line interference
Xpower = [sin(2*pi*50*time);cos(2*pi*50*time);...
    sin(2*pi*60*time);cos(2*pi*60*time)]';
Xnew = Xnew - Xpower*inv(Xpower'*Xpower)*Xpower'*Xnew;

for freq = 1:size(Xtemplate,1)
    Xssvep = squeeze(Xtemplate(freq,:,:));
    % Remove potential SSVEP components at ith frequency
    Xnoise = Xnew - Xssvep*inv(Xssvep'*Xssvep)*Xssvep'*Xnew;
    W = canoncorr(Xnew,squeeze(Xtemplate(freq,:,:)));
    s = Xnew*W; % CCA filtered EEG signal
    sNoise = Xnoise*W; % Noise in filtered signal
    recScore = zeros(Nharmonic,size(s,2));

    recSigma = 1;
    for channel = 1:size(s,2)
        for Nh = 1:Nharmonic
            recP = norm(Xssvep(:,2*(Nh - 1)+1:2*Nh)'*s(:,channel),2)^2;
            % AR(p) for estimate noise power
            if p ~= 0
                [beta,estNoiseVar] = aryule(sNoise(:,1),p); 
                recSigma = (pi*size(Xnew,1)*estNoiseVar/4)/...
                    abs(1+sum(beta(2:end).*exp(-2i*pi*[1:p]*Nh*stimuFreq(freq)/fsample)));
            end
            recScore(Nh,channel) = recP/recSigma;
        end
    end
    score(freq) = mean(mean(recScore));
end

if alpha ~= 0
    scoreNorm = score/sum(score);
    scoreSoftmax = exp(alpha*scoreNorm)/sum(exp(alpha*scoreNorm));
    score = scoreSoftmax;
end
