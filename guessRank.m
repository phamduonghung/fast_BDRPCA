function r = guessRank(M)
[n m] = size(M);
epsilon = nnz(M)/sqrt(m*n);
%S0 = svdecon(M_E);
S0 = svd(M,'eco') ;
%S0 = svds(M_E,100) ;

S1=S0(1:end-1)-S0(2:end);
S1_ = S1./mean(S1(end-10:end));
r1=0;
lam=0.05;
while(r1<=0)
    for idx=1:length(S1_)
        cost(idx) = lam*max(S1_(idx:end)) + idx;
    end
    [v2, i2] = min(cost);
    r1 = max(i2-1);
    lam=lam+0.05;
end

clear cost;
for idx=1:length(S0)-1
    cost(idx) = (S0(idx+1)+sqrt(idx*epsilon)*S0(1)/epsilon)/S0(idx);
end
[v2,i2] = min(cost);
r2 = max(i2);
r = max([r1 r2]);