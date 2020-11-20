clear;
close all;
clc;

% Basic values
lambda_0 = 100;                     % 1/s
r = 45*(10^-9);                     % radius of reciever (m)
d = 5*(10^-7);                    % distance (m)
D = 4.265*(10^-10);                 % Co-efficient of Diffusion (m*m/s)
del_t = 9*(10^-6);                          % Discrete Time length (s)
T = 30*del_t;               % Slot length (s)
L=5;


ri = [0 : 60];
snr_db = 25;
snr_lin = 10^(snr_db/10);
sum_cj = 0;
s_0 = 0;
s_1 = 1;

c_0= snr_lin*2*lambda_0*T;
p_0= (r/d)*(erfc((d-r)/sqrt(4*D*T)));
ntx=c_0/p_0;

for i=1:L
    p(i)=(r/d)*(erfc((d-r)/sqrt(4*D*(i+1)*T))-erfc((d-r)/sqrt(4*D*i*T)));
    c(i)=p(i)*ntx;
    sum_cj = sum_cj + c(i);
end

term1 = (lambda_0*T) + (sum_cj/2);
lam_at_0 = (c_0*s_0) + term1;
lam_at_1 = (c_0*s_1) + term1;
tau_sub_opt = c_0/log(1 + (c_0/term1));
for i=1:length(ri)
    Pe_0(i)= (exp(-lam_at_0) * (lam_at_0)^ri(i)) / factorial(ri(i));
    Pe_1(i)= (exp(-lam_at_1) * (lam_at_1)^ri(i)) / factorial(ri(i));
    if(Pe_0(i)==Pe_1(i))
        tau_opt=ri(i);
    end
end

xyz=figure;
semilogy(ri, Pe_0, 'r' ,ri, Pe_1, 'k')
ylabel('Probability')
xlabel('Received Number of Bits')
legend('Real Dist. for Si=0', 'Real Dist. for Si=1')
title('Emperical Distributions for received Bits (SNR=25dB)')
axis([0 60 10^-6 1])
grid on;
savefig(xyz,'Reproduce-4.fig')