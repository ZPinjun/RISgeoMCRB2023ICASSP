clc; clear; close all

Ptest = -10:5:40;
repeatNum = 80;
data_RMSE = zeros(length(Ptest),4);


%% ================================== mismatch on RIS position
wb = waitbar(0,'Simulating on RIS position mismatch ...');
for i = 1:length(Ptest)
    waitbar(i/length(Ptest),wb);

    P = db2pow(Ptest(i));
    % default setup
    sp = default_setup();
    % update setup
    sp.P = P;
    sp.Pr_error = 0.01*ones(3,1);
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
    Sigma = EFIM_ch^(-1);
    % get pseudo-true parameters and bias
    r_pseudo = get_pseudotrue_CF(sp);
    r_true = [sp.Pu;sp.rho*sp.c];
    bias = norm(r_true(1:3)-r_pseudo(1:3),2)^2;
    % get MCRB and LB
    MCRB = get_MCRB(sp, r_pseudo, Sigma);
    LBM = MCRB + (r_true-r_pseudo)*(r_true-r_pseudo).';
    LB = sqrt(trace(LBM(1:3,1:3)));
    % get RMSE
    [RMSE, ~] = get_RMSE(sp, EFIM_ch, repeatNum);

    data_RMSE(i,1) = PEB;
    data_RMSE(i,2) = RMSE;
    data_RMSE(i,3) = LB;
    data_RMSE(i,4) = sqrt(bias);

end
close(wb);

figure(1)
semilogy(Ptest,data_RMSE(:,1),'k-'); hold on
semilogy(Ptest,data_RMSE(:,2),'bx--'); hold on
semilogy(Ptest,data_RMSE(:,3),'bs-.'); hold on
semilogy(Ptest,data_RMSE(:,4),'b--'); hold off
xlabel('Transmitted power P [dBm]');
ylabel('RMSE');
legend('PEB ','ML-RMSE (pos. mismatch)','LB (pos. mismatch)','Bias (pos. mismatch)');



%% ================================== mismatch on RIS orientation
data_RMSE2 = zeros(length(Ptest),4);

wb = waitbar(0,'Simulating on RIS orientation mismatch ...');
for i = 1:length(Ptest)
    waitbar(i/length(Ptest),wb);

    P = db2pow(Ptest(i));
    % default setup
    sp = default_setup();
    % update setup
    sp.P = P;
    sp.Pr_error = zeros(3,1);
    sp.Or_Euler_error = 0.5*ones(3,1);
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
    Sigma = EFIM_ch^(-1);
    % get pseudo-true parameters and bias
    r_pseudo = get_pseudotrue_CF(sp);
    r_true = [sp.Pu;sp.rho*sp.c];
    bias = norm(r_true(1:3)-r_pseudo(1:3),2)^2;
    % get MCRB and LB
    MCRB = get_MCRB(sp, r_pseudo, Sigma);
    LBM = MCRB + (r_true-r_pseudo)*(r_true-r_pseudo).';
    LB = sqrt(trace(LBM(1:3,1:3)));
    % get RMSE
    [RMSE, ~] = get_RMSE(sp, EFIM_ch, repeatNum);

    data_RMSE2(i,1) = PEB;
    data_RMSE2(i,2) = RMSE;
    data_RMSE2(i,3) = LB;
    data_RMSE2(i,4) = sqrt(bias);

end
close(wb);


figure(1)
semilogy(Ptest,data_RMSE(:,1),'k-'); hold on
semilogy(Ptest,data_RMSE(:,2),'bx--'); hold on
semilogy(Ptest,data_RMSE(:,3),'bs-.'); hold on
semilogy(Ptest,data_RMSE(:,4),'b--'); hold on
semilogy(Ptest,data_RMSE2(:,2),'r+--'); hold on
semilogy(Ptest,data_RMSE2(:,3),'ro-.'); hold on
semilogy(Ptest,data_RMSE2(:,4),'r--'); hold off
xlabel('Transmitted power P [dBm]');
ylabel('RMSE');
legend('PEB ','ML-RMSE (pos. mismatch)','LB (pos. mismatch)','Bias (pos. mismatch)','ML-RMSE (ori. mismatch)','LB (ori. mismatch)','Bias (ori. mismatch)');


