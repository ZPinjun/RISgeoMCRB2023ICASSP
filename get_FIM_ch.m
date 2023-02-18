function FIM_ch = get_FIM_ch(sp)

%% get parameters
K = sp.K;
G = sp.G;
g_d = sp.g_d;
g_r = sp.g_r;
df = sp.df;
c = sp.c;
t_a = sp.t_phi_a;
t_b = sp.t_phi_b;
phi_b_az = sp.phi_b_az;
phi_b_el = sp.phi_b_el;
pAE = sp.pAE_RIS;
N = size(pAE,2);
lambda = sp.lambda;
tau_d = sp.tau_d;
tau_r = sp.tau_r;
sigma = sp.sigma;

%% generate transmitted signal
x = ones(K,1);

%% generate LOS channel
[dtau_d, d_dtau_d] = d(tau_d,K,df,c);
termd = dtau_d .* x;

%% generate NLOS channel
[b_phib, db_phib_dphib_az, db_phib_dphib_el] = b(t_a, t_b, phi_b_az, phi_b_el, pAE, lambda);
[dtau_r, d_dtau_r] = d(tau_r,K,df,c);
rng(1,'simdTwister');
omega = exp(1j*2*pi*rand(N,G));

%% get FIM of channel parameters
FIM_ch = zeros(8,8);
de_tau_d = g_d * d_dtau_d .* x;
de_gd_real = termd;
de_gd_imaginary = 1j*termd;
for i = 1:G
    de_tau_r = g_r * b_phib.' * omega(:,i) * (d_dtau_r .* x);
    de_phib_az = g_r * db_phib_dphib_az.' * omega(:,i) * (dtau_r .* x);
    de_phib_el = g_r * db_phib_dphib_el.' * omega(:,i) * (dtau_r .* x);
    termr = b_phib.' * omega(:,i) * (dtau_r .* x);
    de_gr_real = termr;
    de_gr_imaginary = 1j*termr;
    D = [de_tau_d, de_tau_r, de_phib_az, de_phib_el, de_gd_real, de_gd_imaginary, de_gr_real, de_gr_imaginary];
    FIM_ch = FIM_ch + real(D'*D);
end
FIM_ch = FIM_ch*2/sigma^2;




end

