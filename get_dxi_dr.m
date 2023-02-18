function dxi_dr = get_dxi_dr(pu,pb,pr,Rr)

tub = pu - pb;
tur = pu - pr;
tur_local = Rr.'*tur;
% tau_d
row1 = [(tub/norm(tub,2)).', 1];   % transfer delta_clk from [s] to [m]
% tau_r
row2 = [(tur/norm(tur,2)).', 1];   % transfer delta_clk from [s] to [m]
% phi_b_az
term1 = ( 1 + (tur_local(2)/tur_local(1))^2 )^(-1);
row3 =  [( term1 * ( ( (tur_local(1)*Rr(:,2)) - (tur_local(2)*Rr(:,1)) ) / tur_local(1)^2 ) ).', 0 ];
% phi_b_el
term2 = ( 1 - ( tur_local(3)/norm(tur,2) )^2 )^(-1/2);
row4 = [( term2 * ( Rr(:,3)/norm(tur,2) - tur_local(3)*tur/norm(tur,2)^3 ) ).', 0 ];
dxi_dr = [row1; row2; row3; row4 ];

end

