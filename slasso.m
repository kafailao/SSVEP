function score = slasso(Xnew,Xtemplate,lambda)
% Xnew [samples, channels]
% Xtemplate [freq, samples, channels]
Xsize = size(Xtemplate);
score = zeros(Xsize(1),1);
nXtemplate = reshape(permute(Xtemplate,[3 1 2]),Xsize(1)*Xsize(3),Xsize(2))';
Beta = zeros(size(Xnew,2), Xsize(1)*Xsize(3));
for channel = 1:size(Xnew,2)
     Beta(channel,:) = lasso(nXtemplate,Xnew(:,channel),'lambda',lambda);   
end
% assignin('base','Beta1',Beta);
Beta = reshape(Beta,[size(Xnew,2), Xsize(3), Xsize(1)]);
for freq = 1:Xsize(1)
    score(freq) = sum(sum(abs(squeeze(Beta(:,:,freq)))))/size(Xnew,2);
end
