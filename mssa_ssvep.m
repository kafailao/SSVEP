function score = mssa_ssvep(Xnew,Xtemplate,M)
score = zeros(size(Xtemplate,1),1);
for freq = 1:size(Xtemplate,1)
    [lambda,~,~] = mssa([Xnew,squeeze(Xtemplate(freq,:,:))],M);
    nlabmda = lambda/sum(lambda);
    score(freq) = 1 + sum(bsxfun(@times,nlabmda,log(nlabmda)))/size(Xnew,2);
end
end