function [credibility,data_misfit1,data_misfit2,gamma] = generate_credibility(A,d,m,m_apr,n)

epsilon = 0.001;
d_0 = A*m_apr;
D = A*m;
sigma = var(d);
credibility = zeros(length(m),1);
gamma = 10;
% eta = 1.52;
eta = 0.76;
% eta = 0.792;

if n>20
    eta=eta/exp(n);
end
if eta == 0
    Eta = 1;
else 
    Eta = pinv(eta);
% if n==1
%     eta =0.01;
% end


data_misfit1 = -Eta*norm(D-d_0).^2;
data_misfit2 = eta*norm(D-d,2).^2 ;
% if norm(data_misfit1)<norm(data_misfit2)
%     eta=0;
%     fprintf('wrong\n')
% end
% data_misfit1 = -Eta*norm(D-d_0).^2;
% data_misfit2 = eta*norm(D-d,2).^2 ;
for i = 1:length(m)
%     credibility (i) = 2e+7/(norm(A(:,i)) + epsilon)*exp((data_misfit1+data_misfit2)/(2*sigma^1));
    credibility (i) = 1e+9/(norm(A(:,i)) + epsilon)*exp((data_misfit1+data_misfit2)/(2*sigma^1));
end

end