function r_pseudo = solve_pseudotrue_parameter_ML(r_true, Sigma, pb, pr, Rr, mpr, mRr)

%% =================== use manopt toolbox
% Create the problem structure
manifold_p = euclideanfactory(4,1);
problem.M = manifold_p;
% Define the problem cost function and its Euclidean derivatives
problem.cost = @(r) cost_function(r, r_true, Sigma, pb, pr, Rr, mpr, mRr);
problem.grad = @(r) (manifold_p.egrad2rgrad(r, egrad_function(r, r_true, Sigma, pb, pr, Rr, mpr, mRr)));
% Solve the problem
opt.minstepsize = 1e-30;
opt.maxiter = 80;
%opt.verbosity = 1;
%warning('off', 'manopt:getHessian:approx');
[r_pseudo,~] = trustregions(problem, r_true, opt);


function f_cost = cost_function(r, r_true, Sigma, pb, pr, Rr, mpr, mRr)
    pu = r(1:3,1);
    delta = r(4);
    pu_true = r_true(1:3,1);
    delta_true = r_true(4);
    % get xi & xi_tilde
    xi = get_xi_from_r(pu_true,delta_true,pb,pr,Rr);
    xi_tilde = get_xi_from_r(pu,delta,pb,mpr,mRr);
    f_cost = (xi-xi_tilde).'*Sigma^(-1)*(xi-xi_tilde);
end


function egrad = egrad_function(r, r_true, Sigma, pb, pr, Rr, mpr, mRr)
    pu = r(1:3,1);
    delta = r(4);
    pu_true = r_true(1:3,1);
    delta_true = r_true(4);
    % get xi & xi_tilde
    xi = get_xi_from_r(pu_true,delta_true,pb,pr,Rr);
    xi_tilde = get_xi_from_r(pu,delta,pb,mpr,mRr);
    % df_dxi
    df_dxi = 2*Sigma^(-1)*(xi - xi_tilde);
    % dxi_dr & dxi_tilde_dr
    dxi_tilde_dr = get_dxi_dr(pu,pb,mpr,mRr);

    egrad = - dxi_tilde_dr.'*df_dxi;
end


end

