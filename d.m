function [dtau,d_dtau] = d(tau,K,df,c)

dtau = exp(-1j*2*pi*tau*df*(1:K)./c).';

d_dtau = (-1j*2*pi*df*(1:K).'./c).*dtau;

end

