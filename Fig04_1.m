clc; clear; close all

Ptest_sigma_p = 0.01:0.01:0.06;
repeat = 100;
data_PEB = zeros(length(Ptest_sigma_p),1);
data_LB = zeros(repeat,length(Ptest_sigma_p));

%% ============================================= mismatch on RIS position
% calculate PEB
for i = 1:length(Ptest_sigma_p)
    % default setup
    sp = default_setup();
    % update setup
    sp.Pr_error = Ptest_sigma_p(i)*randn(3,1);
    sp.Or_Euler_error = zeros(3,1);
    sp = update_setup(sp); 
    % get PEB
    FIM_ch = get_FIM_ch(sp);
    EFIM_ch = FIM_ch(1:4,1:4) - FIM_ch(1:4,5:8)*FIM_ch(5:8,5:8)^(-1)*FIM_ch(5:8,1:4);
    Y = get_Jocobian(sp);
    FIM_lo = Y.'*EFIM_ch*Y;
    J = FIM_lo^(-1);
    SPEB = trace(J(1:3,1:3));
    if (SPEB<0)
        SPEB=+inf;
    end
    PEB=sqrt(SPEB);

    data_PEB(i) = PEB;
end
% calculate LB
randerror = randn(3,repeat);
wb = waitbar(0,'Simulating on RIS position mismatch ...');
for j = 1:repeat
    for i = 1:length(Ptest_sigma_p)
        waitbar((j-1)*length(Ptest_sigma_p)/(repeat*length(Ptest_sigma_p)),wb);
        % default setup
        sp = default_setup();
        % update setup
        sp.Pr_error = Ptest_sigma_p(i)*randerror(:,j);
        sp.Or_Euler_error = zeros(3,1);
        sp = update_setup(sp); 
        % get FIM
        FIM_ch = get_FIM_ch(sp);
        EFIM_ch = FIM_ch(1:4,1:4) - FIM_ch(1:4,5:8)*FIM_ch(5:8,5:8)^(-1)*FIM_ch(5:8,1:4);
        Sigma = EFIM_ch^(-1);
        % get true and pseudo-true parameters
        r_pseudo = get_pseudotrue_CF(sp);
        r_true = [sp.Pu;sp.rho*sp.c];
        % get MCRB and LB
        MCRB = get_MCRB(sp, r_pseudo, Sigma);
        LBM = MCRB + (r_true-r_pseudo)*(r_true-r_pseudo).';
        LB = sqrt(trace(LBM(1:3,1:3)));
    
        data_LB(j,i) = LB;
    end
end
close(wb);
data_MCRB_min = min(data_LB);
data_MCRB_max = max(data_LB);
data_MCRB_mean = mean(data_LB);

figure(1)
semilogy(Ptest_sigma_p,data_PEB,'r--', 'LineWidth', 1); hold on
errorbar(Ptest_sigma_p, data_MCRB_mean, data_MCRB_mean-data_MCRB_min, data_MCRB_max-data_MCRB_mean, 'bo-', 'LineWidth', 1); hold off
legend('PEB','LB')
grid on;
xlabel('$\sigma_\mathrm{p}$ [m]','interpreter','latex');
ylabel('RMSE [$m$]','interpreter','latex');


