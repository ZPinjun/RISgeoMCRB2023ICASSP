function [RMSE, Sigma] = get_RMSE(sp, EFIM_ch, repeatNum)

% add noise
xi_true = get_xi_from_r(sp.Pu,sp.rho*sp.c,sp.Pb,sp.Pr,sp.Rr);
[xi, Sigma] = add_noise(xi_true, EFIM_ch, repeatNum);

% get RMSE
data = zeros(1,repeatNum);
pu_true = sp.Pu; 
for i=1:repeatNum
    xi_hat = xi(i,:).';

    r_init = [-2;2;0.2;30]; 
    r_est = position_estimator_ML(r_init, xi_hat, Sigma, sp.Pb, sp.mPr, sp.mRr);
    
    data(i) = norm(pu_true-r_est(1:3,1),2)^2;
end

RMSE = sqrt(mean(data));

end

