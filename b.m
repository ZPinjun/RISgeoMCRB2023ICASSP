function [b_phib, db_phib_dphib_az, db_phib_dphib_el] = b(t_a, t_b, phi_b_az, phi_b_el, pAE, lambda)

% calculate b_phi_b
a_phi_a = exp((1j*2*pi/lambda)*t_a.'*pAE).';
a_phi_b = exp((1j*2*pi/lambda)*t_b.'*pAE).';

b_phib = a_phi_a .* a_phi_b;


% calculate the derivative of b_phi_b
d_az = [-cosd(phi_b_el)*sind(phi_b_az), cosd(phi_b_el)*cosd(phi_b_az), 0];
d_el = [-sind(phi_b_el)*cosd(phi_b_az), -sind(phi_b_el)*sind(phi_b_az), cosd(phi_b_el)];

db_phib_dphib_az = b_phib .* ((-1j*2*pi/lambda) * (d_az*pAE).');
db_phib_dphib_el = b_phib .* ((-1j*2*pi/lambda) * (d_el*pAE).');

end

