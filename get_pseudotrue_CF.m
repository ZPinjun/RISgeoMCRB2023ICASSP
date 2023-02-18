function pseudotrue = get_pseudotrue_CF(sp)

% get parameters
pb = sp.Pb;
Rr = sp.Rr;
mRr = sp.mRr;
pr = sp.Pr;
mpr = sp.mPr;
pu = sp.Pu;
delta = sp.rho*sp.c;

% compute the length of the line segment
a = mRr*Rr.'*(pu-pr);
b = norm(pu-pr,2) + norm(pb-pr,2) - norm(pb-mpr,2) - norm(pb-pu,2);
scale = 0.5*(b^2 - norm(mpr-pb,2)^2) / (a.'*(mpr-pb)+b*norm(a,2));
if scale <0
    scale = nan;
end


% output
pseudotrue = zeros(4,1);
pseudotrue(1:3) = a*scale + mpr;
pseudotrue(4) = delta + norm(pb-pu,2) - norm(pb-pseudotrue(1:3),2);

end

