function sp = default_setup()

%% system setup

% generate BS
sp.Pb = [5, 0, 3].';
% generate UE
sp.Pu = [-2.5 2.5 0].';
% generate RIS
sp.Pr = [0 -5 2.5].';
rng(1,'simdTwister');
sp.Pr_error = 0.03*randn(3,1);
% sp.Pr_error = [0.01,0.01,0.01].';
sp.Or_Euler = [90 0 0].';
sp.Or_Euler_error = 1*randn(3,1);
sp.RIS_dim = [64, 64].';
sp.rho = 1e-7;   % clock offset, in [sec]


%% channel setup

sp.fc = 28e9;   % center frequancy
sp.BW = 400e6;   % Bandwidth
sp.K = 3000;   % # of subcarriers
sp.G = 32;   % # of transmissions
sp.nu = 2;   % path loss component
sp.P = 10;   % average transmit power in [mW]

sp.T0 = 296;                 % Reference temperature (Kelvin)
sp.T = 298.15;               % System temperature (Kelvin), 298.15 Kelvin = 25 celsius
sp.p0 = 1;                   % Standard pressure (atm)
sp.p = sp.p0;                % System pressure (atm)

% constants
sp.c = 2.9979e8;             % Speed of light in vacuum
sp.h = 6.6262e-34;           % Planck constant
sp.Kb = 1.3806e-23;          % Boltzmann constant
sp.R = 8.3144;               % Gas constant
sp.Na = 6.0221e23;           % Avogadro constant
sp.Tstp = 273.15;            % Temperature at standard pressure

end

