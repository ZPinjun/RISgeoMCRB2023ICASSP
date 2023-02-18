function sp = update_setup(sp)

%% get parameters
Pb = sp.Pb;   % position of BS
Pu = sp.Pu;   % position of UE
Pr = sp.Pr;   % position of RIS
Pr_error = sp.Pr_error;
Or_Euler = sp.Or_Euler;   % orientation of RIS
Or_Euler_error = sp.Or_Euler_error;
rho = sp.rho;   % clock offset
c = sp.c;
P = sp.P;
fc = sp.fc;
RIS_dim = sp.RIS_dim;


%% generate RIS rotation matrix
% true orientation
Rr = eul2rotm(reshape(deg2rad(Or_Euler), [1,3]), 'ZYX');
normal_RIS = Rr*[1 0 0].';   % RIS's normal vector in GCS
t_phi_a = Rr.'*(Pb-Pr)/norm(Rr.'*(Pb-Pr),2);
t_phi_b = Rr.'*(Pu-Pr)/norm(Rr.'*(Pu-Pr),2);
% mismatched position & orientation
mPr = Pr + Pr_error;
mRr = eul2rotm(reshape(deg2rad(Or_Euler_error), [1,3]), 'ZYX')*Rr;


%% generate channel parameters: 
% [t_d, t_r, phi_b_az, phi_b_el, g_d, g_r]
dbu = norm(Pb-Pu,2);
dbr = norm(Pb-Pr,2);
dru = norm(Pr-Pu,2);
lambda = c/fc;
tau_d = dbu + rho*c;   % LOS delay, in [m]
tau_r = dbr + dru + rho*c;   % NLOS delay, in [m]
[phi_b_az, phi_b_el] = get_angle_from_dir(t_phi_b);   % AOD from RIS to UE
rng(1,'simdTwister');
g_d = lambda*sqrt(P)*exp(1j*rand(1)*2*pi) / (4*pi*dbu);
g_r = ( (0.5*lambda^2) * sqrt(P) * exp(1j*rand(1)*2*pi) ) / ( (4*pi)^(1.5)*dbr*dru );


%% generate the position of RIS elements
% arrays are deployed on the Y-O-Z plane, x-axis is the normal direction
yrange = (  (1:RIS_dim(1)) - (1+RIS_dim(1))/2  ) * (0.5*lambda);
zrange = (  (1:RIS_dim(2)) - (1+RIS_dim(2))/2  ) * (0.5*lambda);
pAE_RIS = [zeros(1,RIS_dim(1)*RIS_dim(2));
        kron(yrange, ones(1,RIS_dim(1)));
        kron(ones(1,RIS_dim(1)), zrange)];


%% update setup
sp.Rr = Rr;
sp.mPr = mPr;
sp.mRr = mRr;
sp.normal_RIS = normal_RIS;
sp.t_phi_a = t_phi_a;
sp.t_phi_b = t_phi_b;
sp.tau_d = tau_d;
sp.tau_r = tau_r;
sp.phi_b_az = phi_b_az;
sp.phi_b_el = phi_b_el;
sp.g_d = g_d;
sp.g_r = g_r;
sp.lambda = lambda;
sp.pAE_RIS = pAE_RIS;

% generate noise
sp.N0 = sp.Kb*sp.T*1000;   % thermal noise PSD in [mW/Hz]
sp.operationBW = sp.BW;    % Operation bandwidth for Thermal noise
sp.Pn = sp.N0*sp.operationBW;    % thermal noise in [mW]
sp.sigma_in = sqrt(sp.Pn);    % input noise sigma
sp.NoiseFigure = 10;       % noise figure in [dB].
sp.sigma = sqrt(10^(sp.NoiseFigure/10))*sp.sigma_in;

% get frequency fk, k = 1,...,K
df = sp.BW/sp.K;   % Subcarrier bandwidth (Hz)
fstart_sub = sp.fc - sp.BW/2 + df/2;
fstop_sub = fstart_sub + sp.BW - df/2;
sp.fk = fstart_sub:df:fstop_sub;  % Center frequency of each subcarrier (Hz)
sp.df = df;

end

