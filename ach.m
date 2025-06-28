function [x,converse]=ach(P,N0,epsi)
n=N0;
converse=-betaq_low_v2(1-epsi, n, P);
i=1;
taus = linspace(0,1,30).*epsi; 
taus = taus(3:end-2);
for tau=taus
    ka=kappa_inf(tau,P);
   beta1=betaq_up_v2(1-epsi+tau, n, P);
    con(i)=log2(ka)-beta1;
    i=i+1;
end
x= max(con);



