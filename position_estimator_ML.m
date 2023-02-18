function r_est = position_estimator_ML(r_init, xi_hat, Sigma, pb, pr, Rr)

%% =================== use manopt toolbox
% Create the problem structure
manifold_p = euclideanfactory(4,1);
problem.M = manifold_p;
% Define the problem cost function and its Euclidean derivatives
problem.cost = @(r) cost_function(r, xi_hat, Sigma, pb, pr, Rr);
problem.grad = @(r) (manifold_p.egrad2rgrad(r, egrad_function(r, xi_hat, Sigma, pb, pr, Rr)));
% Solve the problem
opt.minstepsize = 1e-30;
opt.maxiter = 100;
% opt.verbosity = 0;
[r_est,~] = trustregions(problem, r_init, opt);




function f_cost = cost_function(r, xi_hat, Sigma, pb, pr, Rr)
    pu = r(1:3,1);
    delta = r(4);
    % get xi
    xi = get_xi_from_r(pu,delta,pb,pr,Rr);
    f_cost = 0.5*(xi_hat-xi).'*Sigma^(-1)*(xi_hat-xi);
end



function egrad = egrad_function(r, xi_hat, Sigma, pb, pr, Rr)
    pu = r(1:3,1);
    delta = r(4);
    % get xi
    xi = get_xi_from_r(pu,delta,pb,pr,Rr);
    % df_dxi
    df_dxi = -Sigma^(-1)*(xi_hat - xi);
    % dxi_dr
    dxi_dr = get_dxi_dr(pu,pb,pr,Rr);
    
    egrad = dxi_dr.'*df_dxi;
end


end

