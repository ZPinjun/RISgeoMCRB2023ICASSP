function r_pseudo = get_pseudotrue_GD(sp, Sigma)

r_true = [sp.Pu;sp.rho*sp.c];

mPr = sp.mPr;
mRr = sp.mRr;

r_pseudo = solve_pseudotrue_parameter_ML(r_true, Sigma, sp.Pb, sp.Pr, sp.Rr, mPr, mRr);

end

