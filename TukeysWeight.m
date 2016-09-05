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