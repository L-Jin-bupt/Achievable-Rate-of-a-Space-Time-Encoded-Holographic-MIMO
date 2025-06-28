clc;
clear all;
close all;
B=10^3; 
f=3*10^9;
c=3*10^8;
lambda=c/f;
L_x=0.5;
L_y=0.2;
epsi=10^-3;
x=10^-2:10^-2:10^-1;
con=1;
SNR=10^3;
for T=x
    %% STE
    N_H=floor(2*pi*T*L_x*L_y*B/(lambda)^2);
    ach_ST(con)=ach(SNR,N_H,epsi)/T;
    disper_H=disper(SNR);
    appro_ST(con)=pi*L_x*L_y*B/(lambda)^2*log2(1+SNR)-sqrt(disper_H*2*pi*L_x*L_y*B/(lambda)^2/T)*qfuncinv(epsi);
    %% TE
    N_T=(2*B*T);
    D_T=(pi*L_x*L_y/(lambda)^2);
    ach_T(con)=D_T*ach(SNR,N_T,epsi)/T;
    appro_T(con)=D_T*(B*log2(1+SNR))-D_T*sqrt(disper_H*2*B/T)*qfuncinv(epsi);
    %% capacity
    capacity(con)=D_T*(B*log2(1+SNR));
    con=con+1;
end

figure
set(gca,'FontName','Times New Rome','FontSize',12);
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultTextFontName','Times New Roman')
plot(x,ach_ST,'-bp',x,appro_ST,'-ro',x,ach_T,'-gv',x,appro_T,'-md',x,capacity,'-k*','LineWidth',1.5);
xlim([10^-2,10^-1]);
xlabel('Transmission duration {\it T} / s');
ylabel('Information rate/ bps')
legend('STE Lem.1','STE Lem.2','TE Lem.3','TE Lem.4','Capacity','Fontsize', 11);
legend('location','southeast');
grid on
hold on