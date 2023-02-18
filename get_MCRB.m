function LB = get_MCRB(sp, r_pseudo, Sigma)

%% get parameters
pb = sp.Pb;
pr = sp.mPr;
Rr = sp.mRr;
pu = r_pseudo(1:3);
delta = r_pseudo(4);
tub = pu - pb;
tur = pu - pr;
tur_local = Rr.'*tur;


%% get matrix A
xi_true = get_xi_from_r(sp.Pu,sp.rho*sp.c,sp.Pb,sp.Pr,sp.Rr);
xi_tilde = get_xi_from_r(pu,delta,pb,pr,Rr);
% ------ the first-order derivative
dxi_dr = get_dxi_dr(pu,pb,pr,Rr);
% ------ the second-order derivative
% tau_d
d2tau_d_pu = eye(3)/norm(tub,2) - tub*tub.'/norm(tub,2)^3;
% tau_r
d2tau_r_pu = eye(3)/norm(tur,2) - tur*tur.'/norm(tur,2)^3;
% phi_b_az
alpha = ( tur_local(1)^2 + tur_local(2)^2 )^(-1);
v = tur_local(1)*Rr(:,2) - tur_local(2)*Rr(:,1);
d_alpha_pu = -( tur_local(1)^2 + tur_local(2)^2 )^(-2) * ( 2*tur_local(1)*Rr(:,1) + 2*tur_local(2)*Rr(:,2) );
d_v_pu = Rr(:,2)*Rr(:,1).' - Rr(:,1)*Rr(:,2).';
d2phi_b_az_pu = v*d_alpha_pu.' + alpha*d_v_pu;
% phi_b_el
alpha = ( norm(tur,2)^2 - tur_local(3)^2 )^(-0.5);
v = Rr(:,3) - tur_local(3)*norm(tur,2)^(-2)*tur;
d_alpha_pu = -0.5*( norm(tur,2)^2 - tur_local(3)^2 )^(-1.5) * ( 2*tur - 2*tur_local(3)*Rr(:,3) );
d_v_pu = tur*( -Rr(:,3)*norm(tur,2)^(-2) + 2*tur_local(3)*norm(tur,2)^(-4)*tur ).' - tur_local(3)*norm(tur,2)^(-2)*eye(3);
d2phi_b_el_pu = v*d_alpha_pu.' + alpha*d_v_pu;
% concatenate them
D2 = zeros(4,4,4);
D2(1:3,1:3,1) = d2tau_d_pu;
D2(1:3,1:3,2) = d2tau_r_pu;
D2(1:3,1:3,3) = d2phi_b_az_pu;
D2(1:3,1:3,4) = d2phi_b_el_pu;
% construct A
A = zeros(4,4);
for i=1:4
    for j=1:4
        d2ij = [D2(i,j,1), D2(i,j,2), D2(i,j,3), D2(i,j,4)].';
        A(i,j) = d2ij.'*Sigma^(-1)*(xi_true-xi_tilde) - dxi_dr(:,i).'*Sigma^(-1)*dxi_dr(:,j);
    end
end


%% get matrix B
B = dxi_dr.'*Sigma^(-1)*( Sigma + (xi_true-xi_tilde)*(xi_true-xi_tilde).' )*Sigma^(-1)*dxi_dr;


LB = A^(-1)*B*A^(-1);


end

