function score = dpls(Xnew, Xtemplate)
% [Input]
% Xnew [samples,channels]: Test EEG data
% Xtemplate [freq,samples,channels]: SSVEP related components, can be
% sinosudial template or EEG template

score = zeros(size(Xtemplate,1),1);
for freq = 1:size(Xtemplate,1)
    Y = squeeze(Xtemplate(freq,:,:));
    W = plsregress(Xnew,Y,size(Y,2));
    B = plsregress(Y,Xnew*W,size(Y,2));
    score(freq) = sum(sum(abs(B)));
end

% Ge, S., Wang, R., Leng, Y., Wang, H., Lin, P., Iramina, K.: A Double-Partial 
% Least-Squares Model for the Detection of Steady-State Visual Evoked Potentials. 
% IEEE Journal of Biomedical and Health Informatics IEEE J. Biomed. Health Inform. 1¡V1 (2016). 


