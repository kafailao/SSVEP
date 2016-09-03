function [p,pvec] = ccaExtend(Xnew, Xtemplate, sinTemplate, mode)
% [Input]
% Xnew [samples,channels]
% Xtemplate [freq,samples,channels]
% sinTemplate [freq,samples,2*Number of Harmonics]
% mode: 'CCA' 'ITCCA' 'Combination3' 'Combination4' 'Combination5'
% [Output]
% p [Number of class,1] vector of correlation coefficient for each class
Wx = zeros(size(Xnew,2),size(sinTemplate,1)); % shape of weight (spatial filter) = [number of channels, number of class]
Wxb = zeros(size(Xnew,2),size(sinTemplate,1));
% Five possible correlation coefficients
p1 = zeros(size(sinTemplate,1),1); % Spatial filter of source template from (template, artifical sin template) apply on new signal and template
p2 = zeros(size(sinTemplate,1),1); % Conventional CCA
p3 = zeros(size(sinTemplate,1),1); % Spatial filter of new signal from (new signal, artifical sin template) apply on new signal and template
%%% Two extra coefficients from PNAS paper
p4 = zeros(size(sinTemplate,1),1); % 
p5 = zeros(size(sinTemplate,1),1);

for i = 1:size(sinTemplate,1)
    if strcmp(mode,'ITCCA')
        Y = squeeze(Xtemplate(i,:,:));
    else
        Y = squeeze(sinTemplate(i,:,:));
    end
        
    % Conventional CCA
    [A,~,R] = canoncorr(Xnew, Y);
    Wx(:,i) = A(:,1);
    p2(i) = R(1);
    if strcmp(mode,'CCA') || strcmp(mode,'ITCCA'), continue; end
    % Two extra coefficients from interSubject paper
    nXtemplate = squeeze(Xtemplate(i,:,:));
    nXtemplate = normalizeSignal(nXtemplate, 1);
    
    temp = corrcoef(Xnew*Wx(:,i),nXtemplate*Wx(:,i));
    p3(i) = temp(2);
    
    [A,~,~] = canoncorr(nXtemplate, Y);
    Wxb(:,i) = A(:,1);
    temp = corrcoef(Xnew*Wxb(:,i),nXtemplate*Wxb(:,i));
    p1(i) = temp(2);
    if strcmp(mode,'Combination3'), continue; end
    %%% Two extra coefficients from PNAS paper
    [A,B,~] = canoncorr(Xnew, nXtemplate);
    temp = corrcoef(Xnew*A(:,1), nXtemplate*A(:,1));
    p4(i) = temp(2);
    if strcmp(mode,'Combination4'), continue; end
    temp = corrcoef(nXtemplate*A(:,1), nXtemplate*B(:,1));
    p5(i) = temp(2);
end
if strcmp(mode,'Combination5')
    p = sign(p1).*p1.^2 + sign(p2).*p2.^2 + sign(p3).*p3.^2 + ...
        sign(p4).*p4.^2 + sign(p5).*p5.^2;
    pvec = [p1 p2 p3 p4 p5];
elseif strcmp(mode,'Combination4')
    p = sign(p1).*p1.^2 + sign(p2).*p2.^2 + sign(p3).*p3.^2 + ...
        sign(p4).*p4.^2;
    pvec = [p1 p2 p3 p4];    
elseif strcmp(mode,'Combination3')
    p = sign(p1).*p1.^2 + sign(p2).*p2.^2 + sign(p3).*p3.^2;
    pvec = [p1 p2 p3];
elseif strcmp(mode,'CCA') || strcmp(mode,'ITCCA')
    p = p2;
    pvec = p;
end

% ['CCA']
% Lin, Z., Zhang, C., Wu, W., Gao, X.: Frequency Recognition Based on Canon-
% ical Correlation Analysis for SSVEP-Based BCIs. IEEE Transactions on Biom-
% edical Engineering IEEE Trans. Biomed. Eng. 53, 2610¡V2614 (2006). 

% ['ITCCA']
% Nakanishi, M., Wang, Y., Wang, Y.-T., Jung, T.-P.: A Comparison Study of 
% Canonical Correlation Analysis Based Methods for Detecting Steady-State 
% Visual Evoked Potentials. PLOS ONE. 10, (2015). 

% ['Combination3']
% Yuan, P., Chen, X., Wang, Y., Gao, X., Gao, S.: Enhancing performances of
% SSVEP-based brain¡Vcomputer interfaces via exploiting inter-subject infor-
% mation. J. Neural Eng. Journal of Neural Engineering. 12, 046006 (2015). 

% ['Combination4']
% Nakanishi, M., Wang, Y., Wang, Y.-T., Mitsukura, Y., Jung, T.-P.: A High-
% Speed Brain Speller Using Steady-State Visual Evoked Potentials. Int. J. 
% Neur. Syst. International Journal of Neural Systems. 24, 1450019 (2014). 

% ['Combination5']
%  Chen, X., Wang, Y., Nakanishi, M., Gao, X., Jung, T.-P., Gao, S.: High-
% speed spelling with a noninvasive brain¡Vcomputer interface. Proceedings 
% of the National Academy of Sciences Proc Natl Acad Sci USA. 112, (2015). 