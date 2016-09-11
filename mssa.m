function [lambda,V,PC] = mssa(X,M)
N = size(X,1);
D = size(X,2);
Y = zeros(N-M+1,D*M);
% create time-delayed embedding of X
X = normalizeSignal(X,1);
for d = 1:D
    for n = 1:N-M+1        
        Y(n,(d-1)*M+1:d*M) = X(n:n+M-1,d);
    end
end

C = Y'*Y/N;
[V,S] = eig(C);
lambda = diag(S);      % extract the diagonal
[lambda,ind]=sort(lambda,'descend'); % sort eigenvalues
V = V(:,ind);             % and eigenvectors

PC = Y*V;