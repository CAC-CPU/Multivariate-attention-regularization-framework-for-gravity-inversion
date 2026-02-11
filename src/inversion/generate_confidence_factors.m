function confFactor = generate_confidence_factors(A,m,d,mapr)
dObs = d;
dApr = A*mapr;
dPre = A*m;
sigma = var(d);
confFactor = zeros(length(m),1);
e = 0.00001;
eta = 0.858;
if eta == 0
    Eta = 1;
else 
    Eta = pinv(eta);

data_misfit = eta*norm(dPre - dObs)^2 - Eta*norm(dPre - dApr)^2;

for i = 1:length(m)
    confFactor(i) = 1/(abs(norm(A(:,i))) + e)*exp(data_misfit/(2*sigma^1));
    if isinf(confFactor(i))
        confFactor(i)=0.001;
    end

end


end