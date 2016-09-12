function [lambda,V,PC,RC] = mssa(X,M)
N = size(X,1);
D = size(X,2);
Y = zeros(N-M+1,D*M);
% create time-delayed embedding of X
X = normalizeSignal(X,1);

for d = 1:D
    pos = (d-1)*M+1:d*M;
    for n = 1:N-M+1        
        Y(n,pos) = X(n:n+M-1,d);
    end
end

C = Y'*Y/N;
[V,S] = eig(C);
lambda = diag(S);      % extract the diagonal
[lambda,ind]=sort(lambda,'descend'); % sort eigenvalues
V = V(:,ind);             % and eigenvectors

remain = ~(cumsum(lambda/sum(lambda)) > 0.9);
if sum(remain) == 0, remain(1) = true; end;
PC = Y*V(:,remain);

RC = zeros(D,N,sum(remain));
for d = 1:D
    for m = 1:sum(remain)
        buf = PC(:,m)*V((d-1)*M+1:d*M,m)'; % invert projection - first channel
        buf = buf(end:-1:1,:);
        for n = 1:N % anti-diagonal averaging
            RC(d,n,m)=mean(diag(buf,-(N-M+1)+n));
        end
    end
end;