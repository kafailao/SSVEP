function score = mssa_ssvep(Xnew,Xtemplate,M)
score = zeros(size(Xtemplate,1),1);
for freq = 1:size(Xtemplate,1)
    % MSI style
    [lambda,~,~] = mssa([Xnew,squeeze(Xtemplate(freq,:,:))],M);
    nlabmda = lambda/sum(lambda);
    score(freq) = 1 + sum(bsxfun(@times,nlabmda,log(nlabmda)))/size(Xnew,2);
    % CCA style
%     [~,~,~,RC] = mssa(Xnew,M);
%     [~,~,r] = canoncorr(squeeze(mean(RC,3))',squeeze(Xtemplate(freq,:,:)));
%     score(freq) = r(1);
end
end