function xi = get_xi_from_r(pu,delta,pb,pr,Rr)

% xi = [tau_d, tau_r, phi_b_az, phi_b_el]
% r = [pu, delta]

tau_d = norm(pb-pu,2) + delta;
tau_r = norm(pb-pr,2) + norm(pr-pu,2) + delta;
t_phi_b = Rr.'*(pu-pr)/norm(Rr.'*(pu-pr),2);
[phi_b_az, phi_b_el] = get_angle_from_dir(t_phi_b);   % AOD from RIS to UE
phi_b_az = deg2rad(phi_b_az);
phi_b_el = deg2rad(phi_b_el);

xi = [tau_d, tau_r, phi_b_az, phi_b_el].';

end

