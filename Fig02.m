clc; clear; close all

%% figure setup
figure(1)
centerX = 300;
centerY = 300;
width =600;
height = 240;
set(gcf,'position',[centerX, centerY,width, height])
FS = 14;


%% ================================================================= beta>0
% default setup
sp = default_setup();
% update setup
sp.Pu = [2 2 0].';
sp.Pr_error = [0.2,0.2,0.2].';
sp.Or_Euler_error = [3,3,3].';
sp = update_setup(sp);
% get FIM
FIM_ch = get_FIM_ch(sp);
EFIM_ch = FIM_ch(1:4,1:4) - FIM_ch(1:4,5:8)*FIM_ch(5:8,5:8)^(-1)*FIM_ch(5:8,1:4);
Sigma = EFIM_ch^(-1);
% get pseudo-true parameters
% ----- the following results should be the same:
r_pseudo_GD = get_pseudotrue_GD(sp, Sigma);   % get pseudo-true parameters by gradient descent method
r_pseudo_CF = get_pseudotrue_CF(sp);   % get pseudo-true parameters by proposed closed-form solution

% plot fig-1
subplot(1,2,1)
mpr = sp.mPr;
pb = sp.Pb;
pu = sp.Pu;
pr = sp.Pr;
beta = norm(pu-pr,2) + norm(pb-pr,2) - norm(pb-mpr,2) - norm(pb-pu,2);
% BS
scatter(pb(1),pb(2),60,'ko','Linewidth', 2);hold on
% RIS
x0 = 0.6;
theta = 0;
aa = [-x0*cosd(theta);-x0*sind(theta)];
bb = -aa;
startpr = pr(1:2) - (bb-aa)/2;
endpr = pr(1:2) + (bb-aa)/2;
plot([startpr(1),endpr(1)], [startpr(2),endpr(2)], 'k','Linewidth', 2); hold on
% mismatched RIS
x0 = 0.6;
theta = 3;
aa = [-x0*cosd(theta);-x0*sind(theta)];
bb = -aa;
startpr = mpr(1:2) - (bb-aa)/2;
endpr = mpr(1:2) + (bb-aa)/2;
plot([startpr(1),endpr(1)], [startpr(2),endpr(2)], 'r','Linewidth', 2);hold on
% true pu
scatter(pu(1),pu(2),200,'gx','Linewidth', 2);hold on
% pseudo-true pu
scatter(r_pseudo_GD(1),r_pseudo_GD(2),200,'rx','Linewidth', 2);hold on
% line
a = sp.mRr*sp.Rr.'*(sp.Pu-sp.Pr);
startp = mpr;
endp = mpr + 5.5*a;
plot([startp(1),endp(1)], [startp(2),endp(2)], 'b--','Linewidth', 0.8);
% heperboloid
x = -10:0.1:10;
y = zeros(1,length(x));
for i = 1:length(x)
    aa = (x(i)-mpr(1))^2 + (r_pseudo_GD(3) - mpr(3))^2;
    bb = (x(i)-pb(1))^2 + (r_pseudo_GD(3) - pb(3))^2;
    y(i) = get_y(aa,bb,beta,pb(2),mpr(2));
end
plot(x,y,'b','Linewidth', 0.8); hold on
% plot edges
plot([-5,5], [-5,-5], 'k-.'); hold on 
plot([5,5], [5,-5], 'k-.'); hold on 
plot([-5,-5], [5,-5], 'k-.'); hold on 
plot([-5,5], [5,5], 'k-.'); hold on 
xlim([-7,7]);
ylim([-7,7]);
legend('$\mathbf{p}_\mathrm{b}$','$\mathbf{p}_\mathrm{r}$', ...
    '$\tilde{\mathbf{p}}_\mathrm{r}$', '$\bar{\mathbf{p}}$',...
    '${\mathbf{p}}_0$','$s_l$', '$s_h$',...
    'interpreter','latex','FontSize',FS);
xlabel('$x$ [m]','interpreter','latex');
ylabel('$y$ [m]','interpreter','latex');
text(-0.3,4,'$\beta>0$','Color','k','FontSize',FS,'HorizontalAlignment','center','interpreter','latex');
box on




%% ================================================================= beta<0
% default setup
sp = default_setup();
% update setup
sp.Pu = [-2 -2 0].';
sp.Pr_error = [0.2,0.2,0.2].';
sp.Or_Euler_error = [3,3,3].';
sp = update_setup(sp);
% get FIM
FIM_ch = get_FIM_ch(sp);
EFIM_ch = FIM_ch(1:4,1:4) - FIM_ch(1:4,5:8)*FIM_ch(5:8,5:8)^(-1)*FIM_ch(5:8,1:4);
Sigma = EFIM_ch^(-1);
% get pseudo-true parameters
r_pseudo_GD = get_pseudotrue_GD(sp, Sigma);   % get pseudo-true parameters by gradient descent method
r_pseudo_CF = get_pseudotrue_CF(sp);   % get pseudo-true parameters by proposed closed-form solution

% plot fig-2
subplot(1,2,2)
mpr = sp.mPr;
pb = sp.Pb;
pu = sp.Pu;
pr = sp.Pr;
beta = norm(pu-pr,2) + norm(pb-pr,2) - norm(pb-mpr,2) - norm(pb-pu,2);
% BS
scatter(pb(1),pb(2),60,'ko','Linewidth', 2);hold on
% RIS
x0 = 0.6;
theta = 0;
aa = [-x0*cosd(theta);-x0*sind(theta)];
bb = -aa;
startpr = pr(1:2) - (bb-aa)/2;
endpr = pr(1:2) + (bb-aa)/2;
plot([startpr(1),endpr(1)], [startpr(2),endpr(2)], 'k','Linewidth', 2); hold on
% mismatched RIS
x0 = 0.6;
theta = 3;
aa = [-x0*cosd(theta);-x0*sind(theta)];
bb = -aa;
startpr = mpr(1:2) - (bb-aa)/2;
endpr = mpr(1:2) + (bb-aa)/2;
plot([startpr(1),endpr(1)], [startpr(2),endpr(2)], 'r','Linewidth', 2);hold on
% true pu
scatter(pu(1),pu(2),200,'gx','Linewidth', 2);hold on
% pseudo-true pu
scatter(r_pseudo_GD(1),r_pseudo_GD(2),200,'rx','Linewidth', 2);hold on
% line
a = sp.mRr*sp.Rr.'*(sp.Pu-sp.Pr);
startp = mpr;
endp = mpr + 5.5*a;
plot([startp(1),endp(1)], [startp(2),endp(2)], 'b--','Linewidth', 0.8);
% heperboloid
x = -10:0.1:10;
y = zeros(1,length(x));
for i = 1:length(x)
    aa = (x(i)-mpr(1))^2 + (r_pseudo_GD(3) - mpr(3))^2;
    bb = (x(i)-pb(1))^2 + (r_pseudo_GD(3) - pb(3))^2;
    y(i) = get_y(aa,bb,beta,pb(2),mpr(2));
end
plot(x,y,'b','Linewidth', 0.8); hold on
% plot edges
plot([-5,5], [-5,-5], 'k-.'); hold on 
plot([5,5], [5,-5], 'k-.'); hold on 
plot([-5,-5], [5,-5], 'k-.'); hold on 
plot([-5,5], [5,5], 'k-.'); hold on  
xlim([-7,7]);
ylim([-7,7]);
% legend('$\mathbf{p}_\mathrm{b}$','$\mathbf{p}_\mathrm{r}$', ...
%     '$\tilde{\mathbf{p}}_\mathrm{r}$', '$\bar{\mathbf{p}}$',...
%     '${\mathbf{p}}_0$','$s_l$', '$s_h$',...
%     'interpreter','latex','FontSize',FS);
xlabel('$x$ [m]','interpreter','latex');
ylabel('$y$ [m]','interpreter','latex');
text(0,4,'$\beta<0$','Color','k','FontSize',FS,'HorizontalAlignment','center','interpreter','latex');
box on
hold off











function yt = get_y(a,b,beta,pb,pr)
range = -10:0.001:10;
vmin = +inf;
ymin = 0;
for j = 1:length(range)
    y = range(j);
    obj = sqrt(a+(y-pr)^2) - sqrt(b+(y-pb)^2) - beta;
    if abs(obj) < vmin
        vmin = abs(obj);
        ymin = y;
    end
end
yt = ymin;
end



