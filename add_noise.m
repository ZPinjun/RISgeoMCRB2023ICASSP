function [measn_all, Sigma] = add_noise(xi_true, FIM, repeatNum)

mu = reshape(xi_true,1,[]);
Sigma = FIM^(-1);
measn_all = mvnrnd(mu,Sigma,repeatNum);

end

